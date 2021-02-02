from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from PIL import Image
from PIL import ImageDraw

# defining the file for the rest of the program using biopython's reader
record = SeqIO.read("Genome.gb", "genbank")

gd_diagram = GenomeDiagram.Diagram(record.id)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

# using the features in the genbank file to create the genome map
for feature in record.features:
    # genbank files have a gene feature that I can separate by checking if it's a gene feature or not
    if feature.type != "gene":
        # this excludes all other features besides genes
        continue
    # checking each feature and separating them into different color sets to make reading the
    # genome map easier
    if len(gd_feature_set) % 2 == 0:
        color = colors.turquoise
    else:
        color = colors.lightgreen
    #setting each feature to the color assigned and adding a label
    gd_feature_set.add_feature(feature, color=color, label=True)

# checking for specific genome sequences and labeling them with their specific labels
for site, name, color in [("GAATTC", "EcoRI", colors.blue), ("GGATCC", "BamHI", colors.red)]:
    index = 0
    while True:
        index = record.seq.find(site, start=index)
        if index == -1: break
        feature = SeqFeature(FeatureLocation(index, index + len(site)))
        gd_feature_set.add_feature(feature, color=color, name=name, label=True, label_size=10, label_color=color)
        index += len(site)

# creating the circular genome map as a png using biopython's built in diagram creation tool
gd_diagram.draw(format="circular", circular=True, pagesize=(20 * cm, 20 * cm),
                start=0, end=len(record), circle_core=0.5)
gd_diagram.write("ToCSV.png", "PNG")

# using Python's Image Library to add a title to the top of the file by opening the png, writing the name of the virus
# and then saving over the file with the title added.
img = Image.open("ToCSV.png")
draw = ImageDraw.Draw(img)
draw.text((200, 10), "Tomato Curly Stunt Virus", (0, 0, 0))
img.save("ToCSV.png")