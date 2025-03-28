import traceback
import json
import boto3
import requests
import os
import logging
from datetime import datetime

ssm_client = boto3.client("ssm")
s3_client = boto3.client("s3", region_name="us-east-2")

# Setup logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def log(message):
    logger.info(message)

class EnvironmentBuild():
    """
    Encapsulated process that allows the environment parameters to be accessed
    by the container
    """
    def __init__(self):
        self.parameter_dictionary = self.__build_paramter_dictionary()

    def __build_paramter_dictionary(self):
        user_cred_esri = os.environ.get("ARCGIS_CREDENTIALS_ARN")
        user_secret_esri = os.environ.get("ARCGIS_SECRET_ARN")
        feature_service_enpoint = os.environ.get("FEATURE_SERVICE_ENDPOINTS")

        user_cred_esri_response = ssm_client.get_parameter(
            Name=user_cred_esri, WithDecryption=True
        )
        user_secret_esri_response = ssm_client.get_parameter(
            Name=user_secret_esri, WithDecryption=True
        )
        feature_service_enpoint_response = ssm_client.get_parameter(
            Name=feature_service_enpoint, WithDecryption=True
        )

        return {
            'client-credential': user_cred_esri_response['Parameter']['Value'],
            'client-secret': user_secret_esri_response['Parameter']['Value'],
            'feature-service-endpoint': feature_service_enpoint_response['Parameter']['Value']
        }

environment_accessor = EnvironmentBuild()

def get_esri_token():
    TOKEN_URL = "https://www.arcgis.com/sharing/rest/oauth2/token"
    client_id = environment_accessor.parameter_dictionary.get('client-credential')
    client_secret = environment_accessor.parameter_dictionary.get('client-secret')

    params = {
        "client_id": client_id,
        "client_secret": client_secret,
        "grant_type": "client_credentials",
        "f": "json"
    }

    response = requests.post(TOKEN_URL, data=params)
    response.raise_for_status()

    token_data = response.json()
    token = token_data.get("access_token")
    if not token:
        raise Exception(f"Token not found in response: {token_data}")
    return token

def update_esri_table(bucket_name=None, s3_event_endpointKeys=None, token=None):
    try:
        log(f"Starting update for key: {s3_event_endpointKeys} in bucket: {bucket_name}")

        feature_service_enpoints = environment_accessor.parameter_dictionary.get('feature-service-endpoint')
        log(type(feature_service_enpoints))
        feature_service_enpoints = json.loads(feature_service_enpoints)
        log(f"Feature service endpoints: {feature_service_enpoints}")

        folder, filename = s3_event_endpointKeys.split('/')
        layer_key = filename.rsplit('.', 1)[0]

        if folder not in feature_service_enpoints:
            raise ValueError(f"Folder '{folder}' not found in feature service mapping.")

        layer_dict = feature_service_enpoints[folder]

        log(f"Layer dict: {layer_dict}")

        s3_response = s3_client.get_object(Bucket=bucket_name, Key=s3_event_endpointKeys)
        table_data = json.loads(s3_response['Body'].read())
        log(f"Table data: {table_data}")


        updated_layers = []

        for service_label, layer_url in layer_dict.items():
            if service_label != layer_key:
                continue

            log(f"Processing {service_label} at {layer_url}")

            delete_resp = requests.post(
                f"{layer_url}/deleteFeatures",
                data={"where": "1=1", "f": "json", "token": token}
            )
            log(f"Delete response: {delete_resp.json()}")
            log(f"Table data: {table_data}")
            features_payload = [{"attributes": row} for row in table_data]
            log(f"Features payload: {features_payload}")
            add_resp = requests.post(
                f"{layer_url}/addFeatures",
                data={
                    "features": json.dumps(features_payload),
                    "f": "json",
                    "token": token
                }
            )
            log(f"Add response: {add_resp.json()}")
            updated_layers.append(service_label)

        if not updated_layers:
            message = f"No matching layers found for key '{layer_key}' in folder '{folder}'."
            log(message)
            return {
                "statusCode": 400,
                "body": message
            }

        success_message = f"Table(s) updated successfully: {updated_layers}"
        log(success_message)
        return {
            "statusCode": 200,
            "body": success_message
        }

    except Exception as e:
        error_trace = traceback.format_exc()
        log(f"Error during table update: {str(e)}\n{error_trace}")
        return {
            "statusCode": 500,
            "error": str(e),
            "traceback": error_trace,
            "input_key": s3_event_endpointKeys
        }

def lambda_handler(event, context):
    results = []
    for record in event['Records']:
        bucket = record['s3']['bucket']['name']
        feature_service_layer_keys = record['s3']['object']['key']
        try:
            log(f"Lambda triggered by key: {feature_service_layer_keys} in bucket: {bucket}")
            token = get_esri_token()
            result = update_esri_table(
                bucket_name=bucket,
                s3_event_endpointKeys=feature_service_layer_keys,
                token=token
            )
        except Exception as e:
            error_trace = traceback.format_exc()
            log(f"Fatal error during processing: {str(e)}\n{error_trace}")
            result = {
                "statusCode": 500,
                "error": f"Failed to process {feature_service_layer_keys}",
                "details": str(e),
                "traceback": error_trace
            }
        results.append(result)

    log(f"Lambda processing complete. Results: {json.dumps(results)}")
    return {
        "statusCode": 207,  # Multi-status to reflect per-record responses
        "results": results
    }
