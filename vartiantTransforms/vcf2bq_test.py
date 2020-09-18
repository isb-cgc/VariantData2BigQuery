import gzip
import gcsfs
import pyarrow
from datetime import datetime 
import csv 
import pandas_gbq 
from tqdm import tqdm as progress
from google.cloud import bigquery
import pandas as pd 
import json
import numpy as np 
import re 
import concurrent.futures

def upload_final_table_version(schema):
    client = bigquery.Client()
    table_id = 'isb-cgc-etl.STAGING.Complete_Target_VCF_rel22'

    # schema = schema
    table = bigquery.Table(table_id,schema=schema)
    table = client.create_table(table)
    print(f"Created table {table.project}, {table.dataset_id}, {table.table_id}")



    job_config = bigquery.QueryJobConfig(
        allow_large_results=True,
        destination=table_id,
    )

    sql = """
        SELECT * 
        FROM `isb-cgc-etl.TARGET_rel22_Vcf.complete_target_vcf`
        """
        
    # Start the query, passing in the extra configuration. 
    query_job = client.query(sql,job_config=job_config)
    query_job.result()
    print(f"Uploaded records to {table.project}, {table.dataset_id}, {table.table_id}")


def upload_to_staging_env(variant_records_file,schema):

    client = bigquery.Client(project='isb-cgc-etl')
    dataset_id = 'TARGET_rel22_Vcf'
    table_id = "test_4"
    full_table_id = "isb-cgc-etl.TARGET_rel22_Vcf.test_4"
    
    table = bigquery.Table(full_table_id,schema=schema)
    table = client.create_table(table)
    print(f"Created table {table.project}, {table.dataset_id}, {table.table_id}")

    dataset_ref = client.dataset(dataset_id)
    table_ref = dataset_ref.table(table_id)
    job_config = bigquery.LoadJobConfig()
    job_config.source_format = bigquery.SourceFormat.CSV 
    job_config.skip_leading_rows = 1

    with open(variant_records_file,'rb') as source_file:
        job = client.load_table_from_file(source_file,table_ref,job_config=job_config)
    
    job.result() # Waits for table load to complete 

    print(f"Loaded {job.output_rows}, rows into {dataset_id}:{table_id}")



def merge_csv_files(file_1,file_2,file_3):
    csv.field_size_limit(10000000)
    csv.field_size_limit()
    with open(file_1,'r') as csv_1, open(file_2,'r') as csv_2, open(file_3,'w') as out_file:
        reader_1 = csv.reader(csv_1)
        reader_2 = csv.reader(csv_2)
        writer = csv.writer(out_file)
        for row_1, row_2 in progress(zip(reader_1,reader_2)):
            writer.writerow(row_1 + row_2)



def create_new_columns(file_1,file_2):
    csv.field_size_limit(10000000)
    csv.field_size_limit()
    with open(file_1) as file_in:
        reader = csv.reader(file_in)
        header = next(reader)
        format_column_index = header.index('FORMAT')
        column_names = set()
        for row in progress(reader):
            cell_information = row[format_column_index]
            column_names.update(cell_information.split(':'))
    column_names = list(column_names)
    num_cols = len(column_names)

    with open(file_1) as file_in:
        with open(file_2,'w') as file_out:
            reader = csv.reader(file_in)
            writer = csv.writer(file_out)
            header = next(reader)
            format_column_index = header.index('FORMAT')
            normal_column_index = header.index('NORMAL')
            tumor_column_index = header.index('TUMOR')
            writer.writerow([f'{name}_Normal' for name in column_names]
                            + [f'{name}_Tumor' for name in column_names])
            for row in progress(reader):
                columns = row[format_column_index].split(':')
                tumor_col_values = row[tumor_column_index].split(':')
                if row[normal_column_index] != '':
                    normal_col_values = row[normal_column_index].split(':')
                    column_indicies = [column_names.index(column) for column in columns]
                    row_out = [''] * (num_cols * 2)
                    for column_index, normal_value, tumor_value in zip(column_indicies,normal_col_values,tumor_col_values):
                        row_out[column_index] = normal_value 
                        row_out[column_index + num_cols] = tumor_value
                    writer.writerow(row_out)
                else:
                    column_indicies = [column_names.index(column) for column in columns]
                    row_out = [''] * (num_cols * 2)
                    for column_index, tumor_value in zip(column_indicies,tumor_col_values):
                        row_out[column_index + num_cols] = tumor_value
                    writer.writerow(row_out)


def generate_dataframe(meta_data, column_headers, records, ref_id,file_url, 
                        project_short_name,file_name,analysis_workflow_type,
                        case_barcode,entity_id,file_1):
    
    """ 
    @paramters 

    Generates one dataframe from VCF file preserving the format of columns 
    and another dataframe containing meta-header, study type, referenge
    genome and sample ID information.

    This function generates 2 dataframes for each VCF file - one dataframe that corresponds
    to all information in the VCF file with the same format of columns and a second datafrarme 
    that contains meta-header information, reference genome version, study type and sample ID
    information.

    @return Pandas Dataframe vcf_df, Pandas Dataframe meta_df 


    """

    # Variant Records Dataframe
    vcf_df = pd.DataFrame(records, columns=column_headers)

    legacy_normal_aliquot_barcode = []
    legacy_tumor_aliquot_barcode = []
    for column in vcf_df.columns:
        if '_aliquot' in column or 'TARGET' in column:
            vcf_df.rename(columns={column:'TUMOR'}, inplace=True)
        elif re.search("[DNA]+_.\d\d$",column) or re.search("T\d$",column):
            legacy_tumor_aliquot_barcode.append(column)
            vcf_df.rename(columns={column:'TUMOR'},inplace=True)          
        elif re.search("N\d$",column):
            legacy_normal_aliquot_barcode.append(column)
            vcf_df.rename(columns={column:'NORMAL'},inplace=True)
        else:
            vcf_df.rename(columns={column:column})
    
    if 'NORMAL' not in vcf_df.columns:
        vcf_df['NORMAL'] = np.nan 
        # vcf_df.columns = column_names
    else:
        column_names = [
        'CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO',
        'FORMAT',
        'TUMOR',
        'NORMAL']     
        vcf_df = vcf_df[column_names]
    
    vcf_df.astype({'POS': 'int64'})
    if len(legacy_tumor_aliquot_barcode) > 0:
        for tumor_id in legacy_tumor_aliquot_barcode:
            vcf_df['legacy_tumor_aliquot_barcode'] = tumor_id
    else:
        vcf_df['legacy_tumor_aliquot_barcode'] = np.nan
    
    if len(legacy_normal_aliquot_barcode) > 0:
        for normal_id in legacy_normal_aliquot_barcode:
            vcf_df['legacy_normal_aliquot_barcode'] = normal_id
    else:
        vcf_df['legacy_normal_aliquot_barcode'] = np.nan

    vcf_df["reference"] = ref_id
    vcf_df["analysis_workflow_type"] = [analysis_workflow_type for wrk_type in range(len(vcf_df))]
    vcf_df["project_short_name"] = [project_short_name for name in range(len(vcf_df))]
    vcf_df["file_gdc_url"] =  [file_url for url in range(len(vcf_df))]
    vcf_df["case_barcode"] = [case_barcode for barcode in range(len(vcf_df))]
    vcf_df["associated_entities__entity_submitter_id"] = [entity_id for sample_id in range(len(vcf_df))]

    return vcf_df 


def parse_vcf(vcf_file):

    """
    @parameters vcf file 
      
    Parses through a VCF file and collects a range of information

    Given each VCF file, this function parses through each row and collects
    meta-header information, Normal and Tumor Sample IDs, reference genome version,
    column headers and records corresponding to each row of information from the file.

    @return meta_data_info[], records[], String ref_id

    """
    
    meta_data_info = []
    records = []
    ref_id = None
    
    for line in vcf_file:
        if line.startswith("##contig"):
            pass
        elif line.startswith("##"):
            meta_data_info.append(line)
        elif line.startswith("#"):
            column_headers = line[1:].strip().split()
        else:
            records.append(line.split())

        if line.startswith("##reference"):
            ref_id = line[-17:].strip()
    
    return meta_data_info, column_headers, records, ref_id

def parse_zipped_vcf(vcf_file):

    """
    @parameters vcf file 
      
    Parses through a VCF file and collects a range of information

    Given each VCF file, this function parses through each row and collects
    meta-header information, Normal and Tumor Sample IDs, reference genome version,
    column headers and records corresponding to each row of information from the file.

    @return meta_data_info[], records[], String ref_id

    """
    
    meta_data_info = []
    records = []
    ref_id = None
    for line in vcf_file:
        if line.decode().startswith("##contig"):
            pass
        elif line.decode().startswith("##"):
            meta_data_info.append(line.decode())
            pass
        elif line.decode().startswith("#"):
            column_headers = line.decode()[1:].strip().split()
        else:
            records.append(line.decode().split())

        if line.decode().startswith("##reference"):
            ref_id = line.decode()[-17:].strip()

    return meta_data_info, column_headers, records, ref_id


def start_process(a_file,project_short_name,file_name, analysis_workflow_type,case_barcode,entity_id,fs,file_1, add_header):  
    if ".gz" in a_file:
        with fs.open(a_file, 'rb') as binary_file: 
            unzipped_file = gzip.GzipFile(fileobj=binary_file)
            meta_data, column_headers, records, ref_id = parse_zipped_vcf(unzipped_file.readlines())
            vcf_df = generate_dataframe(meta_data, column_headers, records, ref_id,a_file,
                            project_short_name,file_name,analysis_workflow_type,case_barcode,entity_id,file_1)        
    else:       
        with fs.open(a_file, 'r') as vcf_file: 
            meta_data, column_headers, records, ref_id = parse_vcf(vcf_file.readlines())
            vcf_df = generate_dataframe(meta_data, column_headers, records, ref_id,a_file,
                            project_short_name,file_name,analysis_workflow_type,case_barcode,entity_id,file_1)
    #print('VCF shape:', vcf_df.shape)
    with open(file_1, 'a') as out_file:
        vcf_df.to_csv(out_file, header=add_header, index=False)


def query_for_table():
    """
    @parameters None 

    SQL query to extract certain columns from bigquery using pandas_gbq.

    This SQL query will select the distinct case_barcodes, file_gdc_id, aliquot_barcode,
    file_gdc_url, analysis_workflow_type, and project_short_name from the...
    rel14_aliquot2caseIDmap, rel14_fileData_current and rel14_GDCfileID_to_GCSurl_NEW
    from ISB-CGC's BigQuery Tables found on Google Cloud. 

    @return pandas_dataframe

    """
    df = pd.io.gbq.read_gbq("""
        SELECT Distinct somatic_hc_snps as file_url, sample_barcode, aliquot_barcode, project_short_name, 
        file_gdc_url, case_barcode,file_gdc_id
        FROM `isb-cgc-etl.TCGA_WGS.tcga_wgs_table`
        LIMIT 2
    """,project_id='isb-cgc-etl',verbose=False,dialect='standard')


    file_urls = list(df["file_url"])
    project_short_name = list(df["project_short_name"])
    aliquot_barcode = list(df["aliquot_barcode"])
    sample_barcode = list(df["sample_barcode"])


    return file_urls,project_short_name,file_name,analysis_workflow_type,case_barcode,entity_id


def schema_with_description():

    with open('~/VariantData2BigQuery/SchemaFiles/tcga_vcf_build_schema_desc.json') as json_file:
        schema = json.load(json_file)

        return schema

def generate_schema_from_csv(file_3):
    
    schema = []
    with open(file_3) as file_in:
        reader = csv.reader(file_in)
        header = next(reader)
        for column_name in header:
            if column_name == 'POS':
                column_dict = {"name":f"{column_name}","type":"integer"}
            else:
                column_dict = {"name":f"{column_name}","type":"string"}
            schema.append(column_dict)

    return schema 

def main():
    startTime = datetime.now()
    file_1 = 'tcga_all.csv'
    file_2 = 'tcga_normal_tumor_columns.csv'
    file_3 = 'tcga_merged.csv'
    bq_schema_with_description = schema_with_description()
    fs = gcsfs.GCSFileSystem(token='google_default')
    max_workers = 70
    
    with open(file_1, 'w') as out_file:
        pass
    
    file_urls,project_short_names,file_names,analysis_workflow_types,case_barcodes,entity_ids = query_for_table()
    
    
    pbar = progress(total=len(file_urls))

    file_urls = iter(file_urls)
    project_short_names = iter(project_short_names)
    file_names = iter(file_names)
    analysis_workflow_types = iter(analysis_workflow_types)
    case_barcodes = iter(case_barcodes)
    entity_ids = iter(entity_ids)


    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        add_header = True
        futures = []
        start_process( 
            next(file_urls), 
            next(project_short_names), 
            next(file_names), 
            next(analysis_workflow_types), 
            next(case_barcodes), 
            next(entity_ids), 
            fs, 
            file_1, 
            add_header)
        pbar.update()
        
        running = set()
        for _, a_file, project_short_name, file_name, analysis_workflow_type,case_barcode,entity_id in zip(
                range(max_workers),
                file_urls, 
                project_short_names, 
                file_names, 
                analysis_workflow_types,
                case_barcodes,
                entity_ids):
            running.add(
                executor.submit( 
                    start_process, 
                    a_file, 
                    project_short_name, 
                    file_name, 
                    analysis_workflow_type, 
                    case_barcode, 
                    entity_id, 
                    fs, 
                    file_1, 
                    False))
        while running:
            done, running = concurrent.futures.wait(running, return_when=concurrent.futures.FIRST_COMPLETED)
            for _ in done:
                pbar.update()
            for _, a_file, project_short_name, file_name, analysis_workflow_type,case_barcode,entity_id in zip(
                    range(len(done)),
                    file_urls, 
                    project_short_names, 
                    file_names, 
                    analysis_workflow_types,
                    case_barcodes,
                    entity_ids):
                running.add(
                    executor.submit( 
                        start_process, 
                        a_file, 
                        project_short_name, 
                        file_name, 
                        analysis_workflow_type, 
                        case_barcode, 
                        entity_id, 
                        fs, 
                        file_1, 
                        False))

    #create_new_columns(file_1,file_2)
    #merge_csv_files(file_1,file_2,file_3)
    
    #schema = generate_schema_from_csv(file_3)
    #upload_to_staging_env(file_3,schema)
    
    upload_final_table_version(bq_schema_with_description)

    print(datetime.now() - startTime)

if __name__ == "__main__":
    main()
