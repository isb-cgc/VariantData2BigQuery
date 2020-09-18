import sys
from json import loads as json_loads
from json import load as json_load
import yaml 
from google.cloud import storage 
from google.cloud import bigquery
from os.path import expanduser
import pandas as pd 
import numpy as np
from alive_progress import alive_bar
from textwrap import dedent
import io
import csv 
import gzip
from google.cloud.exceptions import NotFound
from common_etl.support import confirm_google_vm, parse_genomic_features_file, merge_csv_files,\
                               check_table_existance


def add_labels_and_descriptions(project, full_table_id):
    '''
        @paramaters project, full_table_id

        The function will add in the description, labels, and freindlyName to 
        the published table. 

        @return None 

    '''
    
    client = bigquery.Client(project=project)
    table = client.get_table(full_table_id)

    print('Adding Labels, Description, and Friendly name to table')
    table.description = 'Data was loaded from the ClinVar release created on July 6, 2020. The ClinVar data was loaded from a VCF file downloaded from NCBI and uploaded into BigQuery using the ISB-CGC Variant Converter tool. More details: see https://www.ncbi.nlm.nih.gov/variation/docs/ClinVar_vcf_files/'
    table.friendly_name = 'CLINVAR 20200706 GRCH38'
    assert table.labels == {}
    labels = {"access": "open",
              "data_type":"somatic_mutation",
              "source":"clinvar",
              "reference_genome_0":"hg38",
              "category":"genomic_reference_database",
              "status":"current"}
    table.labels = labels
    table = client.update_table(table, ['description','labels', 'friendlyName'])
    assert table.labels == labels

def publish_table(schema, project, dataset_id, table_id, staging_full_table_id):

    '''
        @parameters schema, project, dataset_id, table_id, staging_full_table_id

        The function will use an SQL query to retrieve the bigquery table from the 
        staging environment and make a push to the production environment including
        the schema description.

        return None

    '''

    client = bigquery.Client(project=project)
    full_table_id = f'{project}.{dataset_id}.{table_id}'

    
    check_table_existance(client,
                          full_table_id,
                          schema)
            
    
    add_labels_and_descriptions(dataset_id,
                                full_table_id)
    

    job_config = bigquery.QueryJobConfig(
        allow_large_results=True,
        destination=full_table_id
    )

    sql = f'''
        SELECT 
            *
        FROM
            {staging_full_table_id}
    '''

    table = client.get_table(full_table_id)
    query_job = client.query(sql, 
                             job_config=job_config)
    query_job.result()
    print(f"Uploaded records to {table.project}, {table.dataset_id}, {table.table_id}")


def upload_to_staging_env(merged_csv, project, dataset_id, table_id):

    client = bigquery.Client(project=project)
    
    df = pd.read_csv(merged_csv, low_memory=False)
    df = df.drop(['INFO'], axis=1)
    table = Table('final_merged.csv', schema='tableschema.json')
    
    
    
    # full_table_id = f'{client.project}.{dataset_id}.{table_id}'
    # job_config = bigquery.LoadJobConfig(write_disposition='WRITE_TRUNCATE')
    # job = client.load_table_from_dataframe(
    #     df, full_table_id, job_config=job_config
    # )
    # job.result()
    
    # table = client.get_table(full_table_id) 
    # print(f'Loaded {table.num_rows} rows and {len(table.schema)} columns to {table_id}')
    



def create_new_columns(file_1, file_2, number_of_lines):

    '''
        @parameters file1, file2, number_of_lines 

        This function will take a csv file with the formated gtf file and 
        parse out the 'attribute' column. The parsed information in the attribute 
        column will be transformed into new columns of their own and be written out
        to a csv file. 

        @return None 

    '''
    header = 1
    number_of_lines_no_header = number_of_lines - header
    
    with open(file_1) as in_file:
        reader = csv.reader(in_file)
        header = next(reader)
        print(header)
        info_column_index = header.index('INFO')
        column_names = []

        with alive_bar(number_of_lines) as bar:
            for row in reader:
                cell_information = row[info_column_index].split(';')
                for cell_info in cell_information:
                    pair = cell_info.split('=')
                    if pair != []:
                        header = pair[0]
                        if header not in column_names:
                            column_names.append(header)
                bar()

        # column_names = list(column_names)
        num_cols = len(column_names)
        print(column_names)
        print(num_cols)

        with open(file_1) as in_file:
            with open(file_2, 'w') as out_file:
                reader = csv.reader(in_file)
                writer = csv.writer(out_file)
                header = next(reader)
                info_column_index = header.index('INFO')
                writer.writerow(column_names)  
                with alive_bar(number_of_lines) as bar:
                    for row in reader:
                        column_indicies = []
                        row_values_list = []
                        cell_information = row[info_column_index].split(';')
                        for cell_info in cell_information:
                                pair = cell_info.split('=')
                                if pair != []:
                                    header = pair[0]
                                    records = pair[1]
                                    column_indicies.append(column_names.index(header))
                                    row_values_list.append(records)
                        row_out = [''] * num_cols
                        for column_index, row_value, in zip(column_indicies,  row_values_list):
                            row_out[column_index] = row_value
                        writer.writerow(row_out)
                        bar()


def parse_vcf(a_file, csv_file):


    with open(csv_file, 'w') as out_file:
        writer = csv.writer(out_file)
        with open('meta_info.txt', 'a') as out_file_2:
            with alive_bar(len(a_file)) as bar:
                for line in a_file:
                    if line.decode().startswith('##'):
                        pass
                    elif line.decode().startswith('##INFO'):
                        out_file_2.write(line.decode())
                    elif line.decode().startswith('#'):
                        column_headers = line.decode()[1:].strip().split()
                        writer.writerow(column_headers)
                    else:
                        line_split = line.decode().strip().split()
                        writer.writerow(line_split)
                    bar()


def count_number_of_lines(a_file):

    '''
        @parameters gff_file

        The number lines in the file will counted and returned
        to keep track of the number of lines being generated for 
        CSV file.

        @return int: number_of_lines
    '''

    number_of_lines = 0 
    with gzip.open(a_file, 'rb') as zipped_file:
        for line in zipped_file:
            if line.decode().startswith("##"):
                pass
            else:
                number_of_lines += 1
    print(f'The number of lines in the file: {number_of_lines}')
    
    return number_of_lines

def schema_with_description(path_to_json):

    '''
        @parameters json_file_path 
        
        Give the file path to the json file to retrieve
        the schema for the table which includes the 
        description as well. 

    '''

    with open(path_to_json) as json_file:
        schema = json_load(json_file)

    return schema

def load_config(yaml_config):
    '''
    ----------------------------------------------------------------------------------------------
    The configuration reader. Parses the YAML configuration into dictionaries
    '''

    yaml_dict = None
    config_stream = io.StringIO(yaml_config)
    try:
        yaml_dict = yaml.load(config_stream, Loader=yaml.FullLoader)
    except yaml.YAMLError as ex:
        print(ex)

    if yaml_dict is None:
        return None, None, None 

    return yaml_dict['files_and_buckets_and_tables'], yaml_dict['steps']



def main(args):
    '''
    Main Control Flow
    Note that the actual steps run are configured in the YAML input! This allows you to
    e.g. skip previously run steps.
    '''

    if not confirm_google_vm():
        print('This job needs to run on a Google Cloud Compute Engine to avoid storage egress charges [EXITING]')
        return


    if len(args) != 2:
        print(" ")
        print(" Usage : {} <configuration_yaml>".format(args[0]))
        return
    
    print("job started")

    with open(args[1], mode='r') as yaml_file:
        params, steps, = load_config(yaml_file.read())


    # variant_call_file_csv = f"~/NextGenETL/files/{params['PARSED_VARIANT_CALL_FILE']}"
    # info_column_split_csv = f"~/NextGenETL/files/{params['INFO_COLUMN_SPLIT_FILE']}"
    # final_merged_csv = f"~/NextGenETL/files/{params['FINAL_MERGED_CSV']}"

    variant_call_file_csv = params['PARSED_VARIANT_CALL_FILE']
    info_column_split_csv = params['INFO_COLUMN_SPLIT_FILE']
    final_merged_csv = params['FINAL_MERGED_CSV']

    # Staging table info for staging env
    staging_project = params['STAGING_PROJECT']
    staging_dataset_id = params['STAGING_DATASET_ID']
    staging_table_id = params['STAGING_TABLE_ID']
    scratch_full_table_id = f'{staging_project}.{staging_dataset_id}.{staging_table_id}'
    
    # Publish table info for production env 
    publish_project = params['PUBLISH_PROJECT']
    publish_dataset_id = params['PUBLISH_DATASET_ID']
    publish_table_id = params['PUBLISH_TABLE_ID']
    path_to_json_schema = params['SCHEMA_WITH_DESCRIPTION']


    if params is None:
        print("Bad YAML load")
        return
    
    schema = schema_with_description(path_to_json_schema)


    if  'count_number_of_lines' in steps:
        print('Counting in the the number of lines in the file')
        number_of_lines = count_number_of_lines(params['FILE_PATH'])

    if 'parse_variant_call_file' in steps:
        print('Processing Variant File')
        with gzip.open(params['FILE_PATH'], 'rb') as zipped_file:
            parse_vcf(zipped_file.readlines(),
                      variant_call_file_csv)

    if 'create_new_columns' in steps:
        print('Creating new columns from the INFO column')
        create_new_columns(variant_call_file_csv, 
                           info_column_split_csv,
                           number_of_lines)

    if 'merge_csv_files' in steps:
        print('Merging CSV files')
        merge_csv_files(variant_call_file_csv,
                        info_column_split_csv,
                        final_merged_csv,
                        number_of_lines)
    
    if 'upload_to_staging_env' in steps:
        print('Uploading table to a staging environment')
        upload_to_staging_env(final_merged_csv,
                              staging_project,
                              staging_dataset_id,
                              staging_table_id)

    if 'publish_table' in steps:
        print('Publishing table')
        publish_table(schema,
                      publish_project,
                      publish_dataset_id,
                      publish_table_id,
                      scratch_full_table_id)




if __name__ == '__main__':
    main(sys.argv)