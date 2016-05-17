for bam in `aws s3 ls s3://amp-alzheimers-rushbroad/Phase1/ | cut -f4 -d' '`; do aws s3 mv s3://amp-alzheimers-rushbroad/Phase1/$bam s3://amp-alzheimers-rushbroad/first_release/$bam; done


for bam in `aws s3 ls s3://amp-alzheimers-rushbroad/Phase2/ --recursive | cut -f4 -d' '`; do aws s3 mv s3://amp-alzheimers-rushbroad/$bam s3://amp-alzheimers-rushbroad/first_release/$bam; done



for fastq in `cut -f4 unmapped.list`; do aws s3 mv s3://amp-mssm-unmapped/$fastq s3://amp-alzheimers-mssm/first_release/unmapped_fastq/$fastq; done