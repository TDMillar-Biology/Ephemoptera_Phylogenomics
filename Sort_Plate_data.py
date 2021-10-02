
## Trevor D. Millar
## This script is used to take in the plate data and make a single csv out file
## File1, File2, Taxa_name (family, genus, species)

import csv
import pdb


with open("./Well3Taxonomy.csv", 'r') as in_file: ## open up our in file
    data = csv.reader(in_file, delimiter = ',') ## read info in in file
            
    rowcount = 0
    out_file = open('Plate_3_key.csv', 'a')
    for row in data:
        if rowcount == 0: ## Skip the first row (header)
            rowcount +=1
            continue
        else:
            ## File_1, File_2, Tax_name
            full_taxonomy = row[4].split('_')
        
            if full_taxonomy[0] == '':
                continue
            else:
                out_file.write('/scratch/kingspeak/serial/u0914152/Ogden_lab/zipped_data/RAPiD-Genomics_F104_' +row[0] + '_i5-502_i7-*_L001_R1_001_val_1.fq.gz' + "," + '/scratch/kingspeak/serial/u0914152/Ogden_lab/zipped_data/RAPiD-Genomics_F104_' + row[0] +'_i5-502_i7-*_L001_R2_001_val_2.fq.gz'  + "," + full_taxonomy[-3] + "_" + full_taxonomy[-2] + "_" + full_taxonomy[-1].replace('.','') + '\n')
                
            
            rowcount += 1



