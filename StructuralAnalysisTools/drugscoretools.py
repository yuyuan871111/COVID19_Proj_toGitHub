import pandas as pd
import os
import re

def drugScoreData2csv(wkdir, FileName, calib_number):
    # Read raw file
    with open(f'{wkdir}/{FileName}.txt', 'r') as F:
        rawtxt_list = F.readlines()
    alascanning_output = []
    for eachline in rawtxt_list:
        m = re.search('^\w+\d+[A-Z]', eachline)
        if not m == None:
            alascanning_output.append(eachline)
    with open(f'{wkdir}/{FileName}_temp.txt', 'w') as F:
        for eachline in alascanning_output:
            F.write(eachline)

    # Read tmp files and give header
    df = pd.read_table(f'{wkdir}/{FileName}_temp.txt', sep = '\t\t', header = None, engine='python')
    mycolnames = ['Residue', 'ddGcalc_kcal/mol', 'DegreeOfPossibleBuriedness', 'Saltbridges']
    df.columns = mycolnames
    # Revise the Residue (remove chain info)
    pd.options.mode.chained_assignment = None
    for index, each in enumerate(df['Residue']):
        m = re.search('(\w\w\w)(\d+)[A-Z]', each)
        residue_calib = int(m.group(2)) + calib_number
        df['Residue'][index] = m.group(1) + str(residue_calib)

    # Write df and remove temp file
    df.to_csv(f'{wkdir}/{FileName}_tidy.csv',index = False)
    os.remove(f'{wkdir}/{FileName}_temp.txt')