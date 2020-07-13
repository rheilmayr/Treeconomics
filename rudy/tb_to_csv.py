import tkFileDialog as fd
import csv
import chardet


def tb_to_csv(tb_file): 
    records = {}

    with open(tb_file, 'rb') as f: 
        reader = csv.reader(x.replace('\0', '') for x in f)
        reader.next()

        for row in reader:  
            # encoding = chardet.detect(''.join(row))['encoding'] 
            # row = [unicode(s.decode(encoding)) for s in row]

            try:           
                key = tuple(row[:7])

                if key not in records: 
                    records[key] = row + [row[-1]]
                else: 
                    records[key][-1] = row[-1]
            except Exception: 
                print e
                print row
                continue
    # print "Records: ", records

    fout_path = fd.asksaveasfilename(title="Save csv as")

    with open(fout_path, 'wb') as fout: 
        writer = csv.writer(fout, delimiter=',')
        
        header = ['Latitude', 'Longitude', 'Elevation', 'Site ID', 'Site Name', 'Species Name', 'Tree ID', 
            'Time Unit', 'Study Start Year BP', 'Study End Year BP', 'Study Start Year', 
            'Study End Year', 'Tree Start Year', 'Tree End Year']

    # Latitude,Longitude,Elevation,Site ID,Site Name,
    # Species Name,Tree ID,Time Unit,Start Year BP,
    # Last Year BP,Start Year,End Year,Tree Year

        writer.writerow(header)
        nprints = 0

        print "\nWriting tb to csv"
        for i, record in enumerate(records): 
            if i == 0: 
                print "\t{} out of {} records".format(i + 1, len(records))

            if nprints >= 1000:                 
                print "\t{} out of {} records".format(i + 1, len(records))
                nprints = 0
            else:
                nprints += 1
            
            try: 
                writer.writerow(records[record])
            except Exception as e: 
                print e
                print records[record]
                continue

        print "\nFinished writing tb to csv"
        




f = fd.askopenfilename(title="Choose tree level records csv")

tb_to_csv(f)
