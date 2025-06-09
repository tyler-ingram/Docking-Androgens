import analyze
import os
import csv

results_dir = "results"
outfile = "results.csv"
dht_center = analyze.get_center('DHT')
with open('results.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['ID', 'Binding Affinity','Tanimoto Similarity', 'Centroid Distance'])
    for filename in os.listdir(results_dir):
        if filename.endswith(".pdbqt"):
            id = filename.split("_")[0]
            tanimoto_similarity = analyze.tanimoto_similarity('DHT', id)
            binding_affinity = analyze.get_binding_affinity(id)
            center = analyze.get_center(id)
            centroid_distance = analyze.calculate_centroid_distance(center, dht_center)
            writer.writerow([id, binding_affinity, tanimoto_similarity, centroid_distance])