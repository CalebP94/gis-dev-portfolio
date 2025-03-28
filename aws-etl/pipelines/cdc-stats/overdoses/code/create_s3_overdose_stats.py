import json
import os
import requests
import boto3
import pandas as pd
import traceback as tb
import io
import logging


TOKEN_URL = "https://www.arcgis.com/sharing/rest/oauth2/token/"
ssm_client = boto3.client("ssm")
lambda_client = boto3.client('lambda')
# Setup logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def log(message):
    logger.info(message)

class EnivronmentAccessor():
    def __init__(self):
        self.environments_dictionary = self.__build_envparam_dictionary()

    def __build_envparam_dictionary(self):
        """
        Build a dictionary of environment credentials and rest api endpoints for use
        in lambda_handlers()
        """
        user_cred_esri = os.environ.get("ARCGIS_CREDENTIALS_ARN")
        user_secret_esri = os.environ.get("ARCGIS_SECRET_ARN")
        scdph_usaboundary_esri = os.environ.get("ARCGIS_USA_POP")
        cdc_covidalloction = os.environ.get("CDC_COVID_ALLOCATION")
        cdc_overdoses = os.environ.get("ARCGIS_CDC_OVERDOSE")
        scdph_countyboundary_esri = os.environ.get("ARCGIS_COUNTIES_POP")
        container_params = os.environ.get("ARCGIS_CONTAINER")
        cdc_apis = os.environ.get("CDC_API_DICT")

        user_cred_esri_response = ssm_client.get_parameter(
            Name=user_cred_esri, WithDecryption=True
        )
        user_secret_esri_response = ssm_client.get_parameter(
            Name=user_secret_esri, WithDecryption=True
        )
        scdph_usaboundary_esri_response = ssm_client.get_parameter(
            Name=scdph_usaboundary_esri, WithDecryption=True
        )
        cdc_covidalloction_response = ssm_client.get_parameter(
            Name=cdc_covidalloction, WithDecryption=True
        )   
        cdc_overdoses_response = ssm_client.get_parameter(
            Name=cdc_overdoses, WithDecryption=True
        )
        scdph_countyboundary_response = ssm_client.get_parameter(
            Name=scdph_countyboundary_esri, WithDecryption=True
        )
        arcgis_container_params = ssm_client.get_parameter(
            Name=container_params, WithDecryption=True
        )
        cdc_apis_response = ssm_client.get_parameter(
            Name=cdc_apis, WithDecryption=True
        )

        return {
            "user_cred_esri": user_cred_esri_response["Parameter"]["Value"],
            "user_secret_esri": user_secret_esri_response["Parameter"]["Value"],
            "scdph_usaboundary_esri": scdph_usaboundary_esri_response["Parameter"]["Value"],
            "scdph_counties_esri": scdph_countyboundary_response["Parameter"]["Value"],
            "cdc_covidalloction": cdc_covidalloction_response["Parameter"]["Value"],
            "cdc_overdoses":cdc_overdoses_response["Parameter"]["Value"],
            "arcgis_container_params":arcgis_container_params["Parameter"]["Value"],
            "cdc_apis": cdc_apis_response["Parameter"]["Value"]
        }
        
environment_accessor = EnivronmentAccessor()

def get_arcgis_credentials():
    """
    Callback to be used in get_esri_token() for access token
    """
    client_id = environment_accessor.environments_dictionary.get("user_cred_esri")
    client_secret = environment_accessor.environments_dictionary.get("user_secret_esri")

    if not client_id:
        raise ValueError("Environment variable 'ARCGIS_CREDENTIALS_ARN' is not set.")
    if not client_secret:
        raise ValueError("Environment variable 'ARCGIS_SECRET_ARN' is not set.")

    return client_id, client_secret


def get_esri_token():
    client_id, client_secret = get_arcgis_credentials()
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

def fetch_esridata_with_token(access_token: str, environment_url:str, build_geom:bool=False):
    if not environment_url and access_token:
        raise ValueError("Environment variable 'ARCGIS_USA_POP' is not set or empty.")
    params = {
        "where": "1=1",
        "token": access_token,
        "outFields": '*',
        "f": "json"
    }
    response = requests.get(environment_url, params=params)
    response.raise_for_status()

    try:
        data = response.json()
        attributes_df = pd.json_normalize([feature["attributes"] for feature in data["features"]])
        attributes_df = attributes_df.drop(columns=[col for col in attributes_df.columns if 'OBJECTID' in col], errors='ignore')

        if build_geom:
            geometry_df = pd.json_normalize([feature.get("geometry", {}) for feature in data["features"]])

            esri_feature_df = pd.concat([attributes_df, geometry_df], axis=1)
        else:
            esri_feature_df = attributes_df

    except json.JSONDecodeError:
        return {"raw_response": response.text}

    return esri_feature_df

def build_cdc_dataframe(url, limit=1000):
    """
    Callback to be used in other funcitons for requesting cdc data
    """
    all_data = []
    offset = 0
    while True:
        paginated_url = f"{url}?$limit={limit}&$offset={offset}"
        response = requests.get(paginated_url)
        response_json = response.json()
        if not response_json:
            break
        all_data.extend(response_json)
        offset += limit
    return pd.DataFrame(all_data)

def normalize_cdc_overdose_data_to_gdb(cdc_api_url:str=None):
    if not cdc_api_url:
        raise ValueError("Environment varibles for cdc apis is not set or empty.")

    cdc_overdose_df = build_cdc_dataframe(cdc_api_url)
    print("Columns:", cdc_overdose_df.columns.tolist())
    print("Data types:\n", cdc_overdose_df.dtypes)

    cdc_overdose_df["provisional_drug_overdose"] = cdc_overdose_df["provisional_drug_overdose"].astype(float)

    cdc_overdose_df["fips_join"] = cdc_overdose_df["fips"].apply(lambda x: '0' + str(x) if len(str(x)) < 5 else str(x))

    aggregated_overdose = cdc_overdose_df.groupby(
        ["year", "countyname", "state_name", "fips_join"]
    )["provisional_drug_overdose"].sum().reset_index()

    return aggregated_overdose

def drop_case_insensitive_duplicates(df):
    seen = {}
    cols_to_drop = []

    for col in df.columns:
        lower_col = col.lower()
        if lower_col in seen:
            cols_to_drop.append(col)
        else:
            seen[lower_col] = col
    log(f"Columns to drop: {cols_to_drop}")
    return df.drop(columns=cols_to_drop)


def standardize_statistics_populationbased(merge_key_df1, merge_key_df2, df1, df2, normalizer, population_column, cdc_sta, statistic_final_column):
    """
    Merge two dataframes based on a common key and perform statistical calculations,
    while retaining all columns from both dataframes.

    Parameters:
    - merge_key_df1 (str): The column name to use as the key for merging in df1.
    - merge_key_df2 (str): The column name to use as the key for merging in df2.
    - df1 (pd.DataFrame): The first dataframe to merge.
    - df2 (pd.DataFrame): The second dataframe to merge.
    - normalizer (int): The value to divide the cdc_sta column by. Standard per x, example 1000.
    - population_column (str): The column name containing the population values.
    - cdc_sta (str): The column name containing the statistic values.
    - statistic_final_column (str): The name of the new column to store the calculated statistic.

    Returns:
    - pd.DataFrame: The merged dataframe with all original columns and the calculated statistic.
    """
    df1 = df1.drop_duplicates(subset=[merge_key_df1])
    df1[merge_key_df1] = df1[merge_key_df1].astype(str)
    df2[merge_key_df2] = df2[merge_key_df2].astype(str)
    
    # Perform the merge, keeping all columns
    merged_df = pd.merge(df1, df2, left_on=merge_key_df1, right_on=merge_key_df2, how='left', suffixes=('_df1', '_df2'))
    # Drop any column that contains 'OBJECTID' in its name
    merged_df.drop(columns=[col for col in merged_df.columns if 'OBJECTID' in col], inplace=True)

    
    # Compute the standardized statistic
    merged_df[statistic_final_column] = (merged_df[cdc_sta] / merged_df[population_column]) * normalizer

    final_df = drop_case_insensitive_duplicates(merged_df)
    
    return final_df


def upload_to_s3(data, file_name):
    try:
        # Ensure input is a DataFrame
        if not isinstance(data, pd.DataFrame):
            raise ValueError("Data must be a pandas DataFrame")

        # Create a byte stream buffer
        buffer = io.BytesIO()

        # Convert DataFrame to JSON and write it to buffer
        json_str = data.to_json(orient="records", indent=4)  # Convert DF to JSON string
        buffer.write(json_str.encode())  # Encode as bytes
        
        # Reset buffer position to start
        buffer.seek(0)

        # Debugging: Print buffer content
        print("File content in buffer:", buffer.getvalue().decode())

        # Initialize S3 client
        s3 = boto3.client("s3", region_name="us-east-2")
        bucket_name = "scdphn-demo-bucket"
        object_key = f"{file_name}"

        try:
            s3.head_object(Bucket=bucket_name, Key=object_key)
            file_exists = True
        except s3.exceptions.ClientError as e:
            if e.response['Error']['Code'] == "404":
                file_exists = False
            else:
                raise

        if file_exists:
            copy_source = {'Bucket': bucket_name, 'Key': object_key}
            
            s3.copy_object(
                CopySource=copy_source,
                Bucket=bucket_name,
                Key=object_key,
                MetadataDirective="REPLACE" 
            )
            print(f"Updated {object_key} in S3 (copied to itself).")
        else:
            s3.put_object(Bucket=bucket_name, Key=object_key, Body=buffer.getvalue())
            print(f"Uploaded {object_key} to S3.")

        return {"status": "Success", "message": f"Uploaded or updated {object_key} in S3"}

    except Exception as e:
        return {"error": str(e)}

def store_json_in_ssm(parameter_name, data_dict):
    event_str = json.dumps(data_dict, indent=2)
    response = ssm_client.put_parameter(
        Name=parameter_name,
        Value=event_str,
        Type="SecureString",
        Overwrite=True
    )
    return response

def lambda_handler(event, context):
    try:
        # Check if the event contains the "ssm_update" key
        if "ssm_update" not in event:
            raise KeyError("Missing 'ssm_update' key in event.")

        param_name = "cdc-container-param"
        resp = store_json_in_ssm(param_name, event)

        ssm_updates = event["ssm_update"]

        if not isinstance(ssm_updates, list):
            raise ValueError("'ssm_update' must be a list of directory-file objects.")

        directory_files_mappings = []
        access_token = get_esri_token()

        cdc_apis_str = environment_accessor.environments_dictionary["cdc_apis"]
        try:
            total_cdcapis = json.loads(cdc_apis_str)
        except json.JSONDecodeError:
            raise ValueError(
                f"CDC_APIS_DICT from SSM is not valid JSON: {cdc_apis_str}"
            )
        stats=[]
        for update in ssm_updates:
            if not isinstance(update, dict):
                raise ValueError("Each item in 'ssm_update' must be a dictionary.")

            # Handle dynamic keys (directory names) in the update object
            for directory, file in update.items():
                if not isinstance(file, str):
                    raise ValueError(f"Files for directory '{directory}' must be a list.")
                if directory in total_cdcapis:

                    cdc_api_url_lookup = total_cdcapis[directory]

                    if "esri_boundary_apiurl" not in event:
                        raise KeyError("Missing 'esri_boundary_apiurl' key in event.")
                    env_dictkey = event["esri_boundary_apiurl"]
                    if not env_dictkey:
                        raise ValueError("No env_dictkey")
                    env_url = environment_accessor.environments_dictionary[env_dictkey]
                    if not env_url:
                        raise ValueError(f"Environment variable '{env_dictkey}' is not set or empty.")

                    esri_feature_data = fetch_esridata_with_token(
                        access_token=access_token,
                        environment_url=env_url
                    )

                    if directory == 'CountiesOverdosePopulationData':
                        cdc_overdoses_df = normalize_cdc_overdose_data_to_gdb(
                            cdc_api_url_lookup
                        )
                        standardized_stats_df = standardize_statistics_populationbased(
                            merge_key_df1="FIPS",
                            merge_key_df2="fips_join",
                            df1=esri_feature_data,
                            df2=cdc_overdoses_df,
                            normalizer=100000,
                            population_column="POPULATION",
                            cdc_sta="provisional_drug_overdose",
                            statistic_final_column="overdose_rate_per_100k"
                        )
                        stats.append({
                                "directory": directory,
                                "file": file,
                                "rowCount": len(standardized_stats_df)
                            })
                        upload_to_s3(standardized_stats_df, f"{directory}/{file}")
                else:
                    raise ValueError(f"Directory '{directory}' not found in CDC_APIS_DICT.")


        return {
            "statusCode": 200,
            "body": stats
        }

    except Exception as e:
        return {
            "statusCode": 500,
            "body": json.dumps({"error": str(e), "traceback": tb.format_exc()})
        }



if __name__ == "__main__":
    lambda_handler(None, None)

