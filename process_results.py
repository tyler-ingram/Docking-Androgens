import analyze
import os
import csv

"""
This script runs all the analysis functions from the analyze.py script on the docking results, and saves the results in a CSV file.
"""

results_dir = "results"
outfile = "results.csv"
dht_center = analyze.get_center('DHT')
with open('results.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    #Headers for the CSV file
    writer.writerow(['ID', 'Binding Affinity','Tanimoto Similarity', 'Centroid Distance'])
    for filename in os.listdir(results_dir):
        if filename.endswith(".pdbqt"):
            #Get the ID from the filename
            id = filename.split("_")[0]
            print(f"Processing {id}...")
            #Run the analysis functions
            tanimoto_similarity = analyze.tanimoto_similarity(id)
            binding_affinity = analyze.get_binding_affinity(id)
            center = analyze.get_center(id)
            centroid_distance = analyze.calculate_centroid_distance(center, dht_center)
            #Write the results to the CSV file
            writer.writerow([id, binding_affinity, tanimoto_similarity, centroid_distance])