# Path to your CSV file
$csvPath = "./developer-m2m_accessKeys.csv"

# Read the CSV and skip the header
$lines = Get-Content $csvPath
$parts = $lines[1] -split ','
Write-Host "Parsed values: $($parts -join ', ')"

$AWS_ACCESS_KEY_ID = $parts[0].Trim()
$AWS_SECRET_ACCESS_KEY = $parts[1].Trim()
$AWS_REGION = $parts[2].Trim()
$AWS_OUTPUT_FORMAT = $parts[3].Trim()
$PROFILE_NAME = "default"

# Create .aws directory if it doesn't exist
$awsDir = "$HOME\.aws"
if (-not (Test-Path $awsDir)) {
    New-Item -Path $awsDir -ItemType Directory | Out-Null
}

# Write credentials file
$credentialsPath = Join-Path $awsDir "credentials"
@"
[$PROFILE_NAME]
aws_access_key_id = $AWS_ACCESS_KEY_ID
aws_secret_access_key = $AWS_SECRET_ACCESS_KEY
"@ | Set-Content -Path $credentialsPath

# Write config file
$configPath = Join-Path $awsDir "config"
@"
[profile $PROFILE_NAME]
region = $AWS_REGION
output = $AWS_OUTPUT_FORMAT
"@ | Set-Content -Path $configPath

Write-Host "AWS CLI configured for profile '$PROFILE_NAME'"

# ECR login
$ecrRegistry = "182399717549.dkr.ecr.us-east-2.amazonaws.com"
Write-Host "Logging into ECR: $ecrRegistry"
aws ecr get-login-password --region $AWS_REGION | docker login --username AWS --password-stdin $ecrRegistry

if ($LASTEXITCODE -eq 0) {
    Write-Host "Docker authenticated with ECR"
} else {
    Write-Host "Docker ECR login failed"
    exit 1
}

docker build -f ./find-near-rehab-center/Dockerfile.dockerfile -t rehabilitation-centers:latest .

aws ecr get-login-password --region us-east-2 | docker login --username AWS --password-stdin 182399717549.dkr.ecr.us-east-2.amazonaws.com

docker tag rehabilitation-centers:latest 182399717549.dkr.ecr.us-east-2.amazonaws.com/rehabilitation-centers:latest

docker push 182399717549.dkr.ecr.us-east-2.amazonaws.com/rehabilitation-centers:latest

# Run the Lambda container locally with credentials passed in
# Write-Host "Running Lambda container locally..."

# docker run `
#   -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID `
#   -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY `
#   -e AWS_DEFAULT_REGION=$AWS_REGION `
#   -p 9000:8080 `
#   $ecrRegistry/arcgis_lambda:latest
