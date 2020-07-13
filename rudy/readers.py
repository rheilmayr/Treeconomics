from core_id_parser import core_id_parser
import re, codecs, chardet, json

CENTI_MM = 0.01
MILLI_MM = 0.001

def get_first(l): 
    if len(l) > 0: 
        return l[0]
    return None

class RwlReader: 
    def __init__(self, package): 
        self.metadata_file = get_first(package['metadata'])
        self.correlation_file = get_first(package['correlation'])
        self.paleodata_file = get_first(package['paleodata'])
        self.site_id = self.site_name = self.species = self.species_id = self.elevation = self.coordinates = self.time_unit = self.year_range = self.year_bp_range = None

        self.units, self.missing_value_id = self.get_units(self.get_content(self.paleodata_file, raw=True))

        self.set_header_from_metadata()
        self.set_header_from_correlation()
        self.set_header_from_rwl()

        if not self.elevation: 
            self.elevation = 0

    def set_header_from_metadata(self): 
        with open(self.metadata_file, 'rb') as mf: 
            metadata = json.load(mf)
            if 'site' in metadata: 
                site_data = get_first(metadata['site'])
                if site_data: 
                    if 'siteName' in site_data: 
                        site_name = site_data['siteName']                 
                        if site_name and len(site_name) and self.site_name == None: 
                            self.site_name = site_name.lower()
                    if 'geo' in site_data: 
                        geo = site_data['geo']
                        if 'geometry' in geo: 
                            geometry = geo['geometry']
                            if 'coordinates' in geometry: 
                                coordinates = geometry['coordinates']
                                if coordinates and self.coordinates == None: 
                                    if len(coordinates) == 2: 
                                        self.coordinates = coordinates
                                    elif len(coordinates) == 4: 
                                        self.coordinates = [coordinates[0], coordinates[2]]
                                        
                        if 'properties' in geo: 
                            properties = geo['properties']
                            if 'maxElevationMeters' in properties: 
                                max_elevation = properties['maxElevationMeters']
                                if max_elevation and len(str(max_elevation)) and self.elevation == None: 
                                    self.elevation = max_elevation
                    if 'paleoData' in site_data: 
                        paleodata = get_first(site_data['paleoData'])

                        if 'timeUnit' in paleodata: 
                            time_unit = paleodata['timeUnit']

                            if time_unit and self.time_unit == None: 
                                self.time_unit = time_unit.strip().lower()

                        if 'earliestYear' in paleodata and 'mostRecentYear' in paleodata: 
                            first_year = paleodata['earliestYear']
                            last_year = paleodata['mostRecentYear']
                            if first_year and last_year and len(str(first_year)) and len(str(last_year)) and self.year_range == None: 
                                self.year_range = [first_year, last_year]
                        
                        if 'earliestYearBP' in paleodata and 'mostRecentYearBP' in paleodata: 
                            first_year_bp = paleodata['earliestYearBP']
                            last_year_bp = paleodata['mostRecentYearBP']
                            if first_year_bp and last_year_bp and len(str(first_year_bp)) and len(str(last_year_bp)) and self.year_bp_range == None: 
                                self.year_bp_range = [first_year_bp, last_year_bp]

                        if 'species' in paleodata: 
                            species_data = get_first(paleodata['species'])

                            if species_data: 
                                if 'speciesCode' in species_data: 
                                    species_code = species_data['speciesCode']
                                    if species_code and len(species_code) and self.species_id == None: 
                                        self.species_id = species_code.lower()
                                if 'commonName' in species_data: 
                                    common_name = species_data['commonName']
                                    if common_name and len(common_name) and self.species == None: 
                                        self.species = ' '.join(common_name).lower()
        
    def set_header_from_correlation(self):         
        content = RwlReader.get_content(self.correlation_file, raw=True) 

        first_year = re.findall(r'Beginning year.*?:.*?([0-9]+)', content)
        last_year = re.findall(r'Ending year.*?:.*?([0-9]+)', content)
        site_name = re.findall(r'Site name.*?:(.*)', content)
        species_info = re.findall(r'Species information.*?:(.*)', content)
        latitude = re.findall(r'Latitude.*?:(.*)', content)
        longitude = re.findall(r'Longitude.*?:(.*)', content)
        elevation = re.findall(r'Elevation.*?:(.*)', content)

        if first_year and last_year: 
            first_year = get_first(first_year).strip()
            last_year = get_first(last_year).strip()
            if not self.year_range:
                self.year_range = [first_year, last_year]
        if site_name: 
            site_name = get_first(site_name).strip().lower()
            if not self.site_name: 
                self.site_name = site_name
        if species_info: 
            species_info = get_first(species_info).split()
            species_id = species_info[0].strip().lower()
            species = ' '.join(species_info[1:]).strip().lower()
            if species_id and species: 
                if not self.species_id: 
                    self.species_id = species_id
                if not self.species: 
                    self.species = species

        if elevation: 
            elevation = get_first(elevation).strip().lower()
            if not self.elevation and elevation.replace('m', ''):
                self.elevation = elevation

    def set_header_from_rwl(self): 
        header = self.get_content(self.paleodata_file, end=3)

        site_id = header[0][:6].strip().lower()
        if site_id and not self.site_id: 
            self.site_id = site_id

        site_name = header[0][9:61].strip().lower()
        if site_name and not self.site_name: 
            self.site_name = site_name

        species_id = header[0][61:65].strip().lower()
        if species_id and not self.species_id: 
            self.species_id = species_id 

        species = header[1][22:40].strip().lower()
        if species and not self.species: 
            self.species = species

        elevation = header[1][40:45].strip().lower().replace('m', '')
        if elevation and not self.elevation: 
            self.elevation = elevation

        coordinates = header[1][47:57].strip().lower()

        year_range = header[1][67:76].strip().split(' ')

        if year_range and len(year_range) > 1 and not self.year_range: 
            self.year_range = year_range

    def get_data(self, test=0): 
        paleodata_rows = self.get_content(self.paleodata_file, start=3)

        def split_row(row): 
            # core_id = row[:6].strip().lower() # row[:8] -> row[:6]
            # decade = row[6:12].strip()  # row[8:12] -> row[6:12]
            start, end = map(int, self.year_range)
            core_id, decade, method = core_id_parser(row[:12], start, end)
            data = row[12:].strip().split()
            
            return core_id, decade, data

        if test: 
            for row in paleodata_rows: 
                yield row
        else:         
            for row in paleodata_rows: 
                try: 
                    core_id, decade, data = split_row(row)
                except: 
                    print(self.paleodata_file, row)
                    break 

                # print core_id, decade, data
                for i, ring_width in enumerate(data):
                    try:  
                        ring_width = int(ring_width)
                    except: 
                        continue

                    try: 
                        decade = int(decade)

                    except: 
                        continue

                    if ring_width != self.missing_value_id:             
                        year = decade + i
                        # print year
                        yield (self.site_id, self.site_name, self.species, self.species_id, self.elevation, 
                                self.coordinates, self.time_unit, self.year_range, self.year_bp_range, core_id, year, 
                                round(ring_width*self.units, 6), decade, row[:12])
            
    @staticmethod
    def get_units(content): 
        if '-9999' in content: 
            return (MILLI_MM, -9999)
        elif '999': 
            return (CENTI_MM, 999)
        else: 
            return 'NA'
    
    @staticmethod
    def detect_encoding(f, nlines=40): 
        i = 0
        with open(f, 'r') as fobj:
            content = []
            for line in fobj: 
                if i > nlines: 
                    break 
                content.append(line)            
                i += 1
        encoding = chardet.detect(''.join(content))['encoding']
        return encoding

    @staticmethod
    def get_content(f, start=0, end=None, raw=False): 
        with codecs.open(f, 'r', encoding=RwlReader.detect_encoding(f)) as fobj: 
            content = fobj.read().replace('\r\n', '\n').replace('\r', '\n')
            
            if raw: 
                return content
            return filter(lambda x: len(x.strip()) > 0, map(lambda x: x.strip(), content.split('\n')[start:end]))