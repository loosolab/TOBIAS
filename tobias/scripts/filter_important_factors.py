#!/usr/bin/env python3
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from argparse import ArgumentParser
import heapq 

#--------------------------------------------------------------------------------------------------------------#
def get_important(args):

    file = args.file_in         #input file bindetect_result.txt
    filter = args.filter        #filter how many binding factors of every condition will be selected
    file_out = args.file_out    #name of output output file
    list_file = []              #contains all lines of file
    new_file = []               #list for the filtered file

    with open(file) as f:   #open bindetect results
        for i in f:
            i = i.strip()
            i = i.split('\t')   #read file tab sapareted
            list_file.append(i)
    
    index_list = [list_file[0].index(i) for i in list_file[0] if '_change' in i]#get the indexs of the columens
    
    important_values = [[max(heapq.nsmallest(filter,[float(a[i]) for a in list_file[1:]])), min(heapq.nlargest(filter,[float(a[i]) for a in list_file[1:]]))]  for i in index_list]

    #important_values contains the maximum and minimum value of the bindingfactor
    for i in list_file[1:]:
        for a,b in zip(index_list, important_values):
            if float(i[a]) >= float(max(b)) or float(i[a]) <= float(min(b)):    #filters if binding value is important
                new_file.append(i)                                              #important lines get append to new list
                print(i[0])                                                     #print stdout for nextflowpipeline
                break                                                           #if line is added for loop jumps to next line 
    
    #build new tab seperater text file 
    book = {i:[] for i in list_file[0]}
    [[book[key].append(value) for key,value in zip(list_file[0], i)] for i in new_file]   #dict for exele wirter key first line value hole line
    df = pd.DataFrame(book)
    df.to_csv(file_out , '\t', index=False)
#--------------------------------------------------------------------------------------------------------#

def main():
    parser = ArgumentParser()
    parser.add_argument("-in", dest="file_in", type=str,   help="Input file")
    parser.add_argument("-filter", dest="filter", type=int,   help="Filter")
    parser.add_argument("-o", dest="file_out", type=str, help="Output file")
    args = parser.parse_args()

    get_important(args)

#--------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    main()









