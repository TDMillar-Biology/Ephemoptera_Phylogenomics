## Trevor D. Millar
## Im writting this script to iterate through my output directories and find which spades assemblies failed. 
## The spades summary report output tells you which assemblies were successful.  


import os
failures = []
no_file = []
successes = []
to_check = []

working_directory = '/uufs/chpc.utah.edu/common/home/u0914152/Spades' ## This is where all of my spades output directories are kept
with open('Plate_3_key.csv', 'r') as in_file:
    all_lines = in_file.readlines()

for line in all_lines:
    to_check.append(line.split(',')[2].strip())



for species in to_check:
    try:
        
        
        check_for = species + '_targetsFULL_ORTHO.fasta'
        
        if check_for in os.listdir(species):
            successes.append(species)
        else:
            failures.append(species)
    except:
        no_file.append(species)

with open('Spades_Summary_report.txt', 'w') as out_file:
    for species in successes:
        out_file.write(species + ', successfull \n')
    out_file.write('#######################################\n')
    for item in failures:
        out_file.write(item + ', failed \n')
    out_file.write('#######################################\n')
    for item in no_file:
        array_num = (to_check.index(item) + 1)
        out_file.write(item + ', not found array = ' + str(array_num) +'\n')



    

