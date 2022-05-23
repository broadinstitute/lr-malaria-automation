import os
import re
import hashlib
import math

import pandas as pd
import firecloud.api as fapi
import numpy as np

#from google.cloud import bigquery
from google.cloud import storage
from google.api_core.exceptions import NotFound

from collections import OrderedDict

import xmltodict
import pprint

import argparse


def load_summaries(gcs_buckets):
    storage_client = storage.Client()
    schemas = OrderedDict()

    ts = {}
    for gcs_bucket in gcs_buckets:
        blobs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket), prefix="inputs")

        for blob in blobs:
            if blob.name.endswith(".fast5") and "fail" not in blob.name:
                gcs_path = os.path.dirname(os.path.dirname(blob.name))
                
                ts[gcs_path] = blob.time_created

    return ts


def upload_sample_set(namespace, workspace, tbl):
    # delete old sample set
    ss_old = fapi.get_entities(namespace, workspace, f'sample_set').json()
    sample_sets = list(map(lambda e: e['name'], ss_old))
    f = [fapi.delete_sample_set(namespace, workspace, sample_set_index) for sample_set_index in sample_sets]

    # upload new sample set
    ss = tbl.filter(['participant'], axis=1).drop_duplicates()
    ss.columns = [f'entity:sample_set_id']
    
    b = fapi.upload_entities(namespace, workspace, entity_data=ss.to_csv(index=False, sep="\t"), model='flexible')
    if b.status_code == 200:
        print(f'Uploaded {len(ss)} sample sets successfully.')
    else:
        print(b.json())
    
    # upload membership set
    ms = tbl.filter(['participant', 'entity:sample_id'], axis=1).drop_duplicates()
    ms.columns = [f'membership:sample_set_id', f'sample']
    
    c = fapi.upload_entities(namespace, workspace, entity_data=ms.to_csv(index=False, sep="\t"), model='flexible')
    if c.status_code == 200:
        print(f'Uploaded {len(ms)} sample set members successfully.')
    else:
        print(c.json())


def merge_tables(tbl_old, tbl_new):
    if tbl_old is not None:
        outer_tbl = pd.merge(tbl_old, tbl_new, how='outer', sort=True, indicator=True)
    else:
        outer_tbl = tbl_new

    hs = []
    for l in list(outer_tbl['entity:sample_id'].unique()):
        g = outer_tbl.loc[outer_tbl['entity:sample_id'] == l].sort_values('_merge')

        if len(g) == 1:
            hs.append(g.iloc[0].to_dict())
        else:
            print(g)
            h = {}
            for col_name in list(outer_tbl.columns):
                q = g[col_name]
                v = q.where((q != 'None') & (q != 'nan')).dropna()
                h[col_name] = v.iloc[0] if len(v) > 0 else ''

            hs.append(h)

    joined_tbl = pd.DataFrame(hs)

    if '_merge' in joined_tbl:
        del joined_tbl['_merge']
    c = list(joined_tbl.columns)
    c.remove("entity:sample_id")
    c = ["entity:sample_id"] + c
    joined_tbl = joined_tbl[c]

    return joined_tbl


def main():
    parser = argparse.ArgumentParser(description='Update Terra workspace sample table', prog='update_nanopore_tables')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-r', '--run', action='store_true', help="Turn off the default dry-run mode")
    parser.add_argument('buckets', metavar='B', type=str, nargs='+', help='GCS buckets to scan')
    args = parser.parse_args()

    print(fapi.whoami())

    print("test1")
    ent_old = fapi.get_entities(args.namespace, args.workspace, 'sample').json()
    print("test2")

    if len(ent_old) > 0:
        tbl_old = pd.DataFrame(list(map(lambda e: e['attributes'], ent_old)))
        tbl_old["entity:sample_id"] = list(map(lambda f: f['name'], ent_old))

    configs = {
        'FLO-MIN106': {'SQK-LSK109': 'dna_r9.4.1_450bps_sup.cfg'},
        'FLO-MIN111': {'SQK-LSK109': 'dna_r9.4.1_450bps_sup.cfg',
                       'SQK-LSK110': 'dna_r10.3_450bps_sup.cfg'},
        'FLO-MIN112': {'SQK-LSK112': 'dna_r10.4_e8.1_sup.cfg'},
        'FLO-PRO002': {'SQK-LSK109': 'dna_r9.4.1_450bps_sup_prom.cfg'}
    }

    ts = load_summaries([default_bucket])

    storage_client = storage.Client()

    tbl_rows = []
    for k,v in ts.items():
        if k != "inputs":
            f = default_bucket + "/" + k

            fs, ss = "unknown", "unknown"
            d = {
                'instrument': 'unknown', #instrument=MC-110675
                'position': 'unknown', #position=MC-110675_0
                'protocol_group_id': 'unknown', #protocol_group_id=coi2_17may2021
                'flow_cell_id': 'unknown', #flow_cell_id=FAO99587
                'sample_id': 'unknown', #=no_sample
                'protocol': 'unknown', #=sequencing/sequencing_MIN106_DNA:FLO-MIN106:SQK-LSK109
                'basecalling_enabled': 'unknown', #==1
            }

            blobs = storage_client.list_blobs(re.sub("^gs://", "", default_bucket), prefix=k)
            for b in blobs:
                if re.search("\/final_summary.*.txt", b.name):
                    fs = default_bucket + "/" + b.name

                    q = b.download_as_text().splitlines()

                    for a in q:
                        if '=' in a:
                            k2,v2 = a.split("=")

                            if v2 != "":
                                d[k2] = v2

                if re.search("\/sequencing_summary.*.txt", b.name):
                    ss = default_bucket + "/" + b.name

            basecalling_model = ''
            a = re.sub('sequencing/sequencing_', '', d['protocol']).split(":")
            if len(a) == 3:
                inst_type, flowcell_type, kit_type = a

                basecalling_model = configs[flowcell_type][kit_type]
                    
            tbl_rows.append([
                hashlib.md5(f.encode("utf-8")).hexdigest(),
                f,
                fs,
                ss,
                d['instrument'],
                d['position'],
                d['protocol_group_id'],
                d['flow_cell_id'],
                d['protocol'],
                str(v),
                d['sample_id'],
                basecalling_model,
                "NA",
                "none",
                os.path.basename(f)
            ])
        
    tbl_header = ["entity:sample_id", "fast5_dir", "final_summary", "sequencing_summary", "instrument", "position",
                  "protocol_group_id", "flow_cell_id", "protocol", "upload_date", "sample_name", "basecalling_model",
                  "barcode_kit", "notes", "fast5_dir_basename"]
        
    tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)
    tbl_new = tbl_new.sort_values(['fast5_dir_basename', 'upload_date'])
    tbl_new = tbl_new.groupby('fast5_dir_basename').first().reset_index().drop(columns='fast5_dir_basename')

    outer_tbl = merge_tables(tbl_old, tbl_new)
    outer_tbl = outer_tbl.sort_values(['final_summary', 'flow_cell_id', 'reference']).groupby('flow_cell_id').first().reset_index()

    final_tbl = outer_tbl
    columns_reordered = ['entity:sample_id'] + list(filter(lambda x: x != 'entity:sample_id', list(final_tbl.columns)))
    final_tbl = final_tbl[columns_reordered]

    pd.set_option('max_columns', 200)
    pd.set_option("max_colwidth", None)

    print(final_tbl)

    #a = fapi.upload_entities(namespace, workspace, entity_data=final_tbl.to_csv(index=False, sep="\t"), model='flexible')

    #if a.status_code == 200:
        #print(f'Uploaded {len(final_tbl)} rows successfully.')
    #else:
        #print(a.json())

    # old
    #s = update_sample_table(args.namespace, args.workspace, args.buckets, args.project)
    #ss, nms = update_sample_set_table(args.namespace, args.workspace, s)

    #if args.run:
        #upload_tables(args.namespace, args.workspace, s, ss, nms)


if __name__ == "__main__":
    main()
