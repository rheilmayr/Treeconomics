from finders import rwl_finder
# from readers import RwlReader

from tkinter import filedialog as fd
import csv, time

def core_id_parser_bf(s, start, end):  
    s_reverse = s[::-1].strip()
    year = []

    for char in s_reverse: 
        if not char.isdigit(): 
            if char == '-':
                pass
            else: 
                break
        
        if len(year) >= len(str(start)): 
            try:
                int(''.join([char] + year))
            except: 
                break 

            if start < int(''.join([char] + year)) < end: 
                year.insert(0, char)
            else: 
                break
        else: 
            year.insert(0, char)

    year = ''.join(year)
    core_id = s[:s.index(year)]
    return core_id, year, 5

def core_id_parser(row, start, end): 
    if row[7] == "." or row[7] == "0" or row[7].isalpha() or row[7] == " ":         
        return row[:8], row[8:].strip(), 0
    elif row[8] == " " and row[8: ].strip().isdigit():
        return row[:8], row[8:].strip(), 1
    elif row[7] == "'":
        year = row.split("'")[-1]
        core_id = row[:row.index("'{}".format(year)) + 1].strip()
        return core_id, year, 2
    
    elif row[7] == "-" and (row[6] == " " or row[6].isalpha()): 
        year = "-{}".format(row.split("-")[-1])
        core_id = row[:row.index(year)].strip()
        return core_id, year, 3
    elif row[7].isdigit() and row[7:].isdigit(): 
        partial_row = row[7:]
        if start < int(partial_row) < end: 
            return row[:7], partial_row, 4
        else:
            return row[:8], partial_row[1:], 4
    else: 
        return core_id_parser_bf(row, start, end)
    
# NEED TO DEAL WITH boi508a12000

# l = [
#     " 9_Col  1800"
# ]


# for i in l: 
#     print core_id_parser(i, 476, 1999)




    # deal with TPDS03C21750 
    # deal with 112125911754
    # deal with 11211174 980
    # deal with boi206b11040
    # deal with  776052 1920
    # deal with SF 31 E'2000
    # deal with ggeo0501b184
    # deal with  788031 1970
    # deal with DLH232B1 990
    # deal with  1_Col  1700
    # deal with  PLK18A	1990
    # deal with  15502  1919
    # deal with dj2--L061660
    # deal with GKMts_121990
    # deal with MNP312M-1370
    # deal with rlztl-092000
    # deal with shut26-21711
    # deal with TC80_402 200
    # deal with  L1V05  1970
    # deal with a9_09me_ 682
    # deal with MW321x -2420
    # deal with PP356a -1120
    # deal with 153  1960  1
    # deal with G024   -2190
    # deal with TPX1@bs'1510
    # deal with BMCH07-21910
    # deal with UC06-11817 6
    # deal with  2 281  1950
    # deal with 72-10-161937















# rwl_dir = fd.askdirectory(title="Choose directory with RWLs")
# # rwl_dir = r"E:\gdp_temperature_project\results\treering_data_width"

# # with open("core_ids_normal.csv", "wb") as core_ids_normal, open("core_ids_weird.csv", "wb") as core_ids_weird:
# with open("core_id_parser_{}.csv".format(time.time()), "wb") as out_csv:
#     writer = csv.writer(out_csv, delimiter=",")
#     # writer_normal = csv.writer(core_ids_normal, delimiter=",")
#     # writer_weird = csv.writer(core_ids_weird, delimiter=",")
#     packages = [i for i in rwl_finder(rwl_dir)]
#     for i, package in enumerate(packages): 
#         print "{} out of {}".format(i + 1, len(packages))
#         package_paleodata_files = package['paleodata']

#         for paleodata_file in package_paleodata_files: 
#             package_copy = {
#                 'correlation': package['correlation'], 
#                 'paleodata': [paleodata_file],
#                 'metadata': package['metadata']
#             }

#             try: 
#                 reader = RwlReader(package_copy)
#             except Exception as e: 
#                 print "RWL reader error", e, paleodata_file
#                 # writer.writerow([paleodata_file, e])
#                 continue

#             for row in reader.get_data(test=1): 
#                 try: 
#                     row = row[:12]
#                     start, end = map(int, reader.year_range)

#                     core_id, year, method = core_id_parser(row, start, end)

#                     writer.writerow([paleodata_file, row, core_id, year, method])

#                 except Exception as e: 
#                     print "core_id_parser error", e
