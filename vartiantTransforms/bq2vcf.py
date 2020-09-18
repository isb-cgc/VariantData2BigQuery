from google.cloud import bigquery
import sys
from oauth2client import client
import pandas as pd
import gcsfs 
import gzip

def create_vcf_from_bigquery(bq_result, column_names, meta_header_info,file_name):
    bq_result =  bq_result.drop(bq_result.columns[-1],axis=1)
    with open(file_name,'w') as vcf_file: 
        for meta_header in meta_header_info:
            vcf_file.write(f'{meta_header}')     
        vcf_file.write('#')
        for header in column_names:
            vcf_file.write(f'{header}\t')
        vcf_file.write('\n')
        
        bq_result.to_csv(vcf_file,sep='\t',index=False,header=False)
    

def get_bigquery_schema(vcf_column_names):

    client = bigquery.Client(project='FILL THIS IN')
    table_id = 'FILL THIS IN'
    table = client.get_table(table_id)  
    schema = table.schema

    meta_header_info = []
    for row in schema:
        if row.name in vcf_column_names:
            meta_header_info.append((f'##{row.name}=<Description="{row.description}">'))

    return meta_header_info

def get_meta_header_from_url(file_url):
    fs = gcsfs.GCSFileSystem(project='FILL THIS IN')
    meta_data_info = []
    with fs.open(file_url[0],'rb') as vcf_file:
        unzipped_file = gzip.GzipFile(fileobj=vcf_file)
        for line in unzipped_file:
            if line.decode().startswith('##contig'):
                pass
            elif line.decode().startswith('##'):
                meta_data_info.append(line.decode())
    
    return meta_data_info

def get_bigquery_table():
    client = bigquery.Client(project='FILL THIS IN')
    query_job = ('''
        SELECT 
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,NORMAL,TUMOR,file_gdc_url as url
        FROM 
            `isb-cgc-etl.TCGA_WGS.somatic_hc_snps_vcf`
        WHERE 
            POS BETWEEN 30000 AND 400000
            AND FILTER = 'PASS'
        ORDER BY
            CHROM
        LIMIT 500
    ''')
    
    dataframe = (
        client.query(query_job).to_dataframe()
    )
    # print(results.schema)
 
        
    return dataframe 


def vcf_column_type(user_input):
    if user_input.lower() == 'controlled':
        fixed_vcf_columns_controlled = [
            'CHROM',
            'POS',
            'ID',
            'REF',
            'ALT',
            'QUAL',
            'FILTER',
            'INFO',
            'FORMAT',
            'NORMAL',
            'TUMOR'
        ]
        return fixed_vcf_columns_controlled

    else:
        fixed_vcf_column_names_open = [
            'CHROM',
            'POS',
            'ID',
            'REF',
            'ALT',
            'QUAL',
            'FILTER',
            'INFO',
            'FORMAT',
        ]
        return fixed_vcf_column_names_open


def main():
    
    user_input_for_file_type = sys.argv[1] # Controlled or Open Acces
    results = get_bigquery_table()
    meta_header_info = get_meta_header_from_url(results['url'].unique())
    vcf_column_names = vcf_column_type(user_input_for_file_type)
    # meta_header_info = get_bigquery_schema(vcf_column_names)
    create_vcf_from_bigquery(results,vcf_column_names,meta_header_info,'test_vcf.vcf')


if __name__ == '__main__':
    main()