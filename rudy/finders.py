import os

def rwl_finder(path): 
    directories = (os.path.join(path, d) for d in os.listdir(path))

    for directory in directories: 
        package = {'paleodata' : [], 'metadata' : [], 'correlation' : []}

        for root, dirs, files in os.walk(directory): 
            if root.endswith('metadata'): 
                package['metadata'].extend(map(lambda p: os.path.join(root, p), files))
            if root.endswith('correlation-stats'): 
                package['correlation'].extend(map(lambda p: os.path.join(root, p), files))

            if 'measurements' in root and all([f.endswith('.rwl') for f in files]):         
                package['paleodata'].extend(map(lambda p: os.path.join(root, p), files))
        
        if len(package['paleodata']) > 0: 
            yield package
        else: 
            print("ERROR ", directory, package)
            continue

# path = fd.askdirectory(title="Choose directory with RWL files")

# for package in rwl_finder(path): 
#     print package