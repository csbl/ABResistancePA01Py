import GEOparse

gse = GEOparse.get_GEO(filepath='./GSE30021_family.soft')

print gse.gsms['GSM743001']