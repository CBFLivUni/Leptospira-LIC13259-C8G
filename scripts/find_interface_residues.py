import gemmi
import os
import csv
import argparse
from collections import Counter

# Function to get residue details
def get_residue_details(entity):
    return f"{entity.residue.name} {entity.residue.seqid.num} {entity.chain.name}"

# Function to find contacts between specified chains in a PDB file
def find_contacts_in_pdb(pdb_file, chain1_id, chain2_id):
    # Load the structure
    structure = gemmi.read_structure(pdb_file)

    # Perform neighbor search on the first model
    ns = gemmi.NeighborSearch(structure[0], structure.cell, 5).populate()

    # Set up contact search with a cutoff distance of 4.0 Å
    cs = gemmi.ContactSearch(4.0)

    # Find contacts
    results = cs.find_contacts(ns)

    # Set to store unique contacts
    unique_contacts = set()

    # Loop through results and add unique residue details for each contact between specified chains
    for contact in results:
        chain1 = contact.partner1.chain.name
        chain2 = contact.partner2.chain.name
        if (chain1 == chain1_id and chain2 == chain2_id) or (chain1 == chain2_id and chain2 == chain1_id):
            partner1 = get_residue_details(contact.partner1)
            partner2 = get_residue_details(contact.partner2)
            # Add tuple of residue details to the set
            unique_contacts.add((partner1, partner2))

    return unique_contacts

def main(pdb_directory, chain1_id, chain2_id):
    # Output directory for .csv files inside the specified pdb_directory
    output_directory = os.path.join(pdb_directory, "contacts_csv")

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Counter to aggregate contact counts
    contact_counter = Counter()

    # Process each .pdb file in the directory
    for filename in os.listdir(pdb_directory):
        if filename.endswith(".pdb"):
            pdb_file = os.path.join(pdb_directory, filename)
            # Find contacts in the PDB file
            contacts = find_contacts_in_pdb(pdb_file, chain1_id, chain2_id)
            # Write contacts to a .csv file
            csv_file = os.path.join(output_directory, f"{os.path.splitext(filename)[0]}_contacts.csv")
            with open(csv_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Residue 1', 'Residue 2'])
                for contact in contacts:
                    writer.writerow([contact[0], contact[1]])
                    # Update the contact counter
                    contact_counter[(contact[0], contact[1])] += 1

    # Write the summary .csv file
    summary_csv_file = os.path.join(output_directory, "summary_contacts.csv")
    with open(summary_csv_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Residue 1', 'Residue 2', 'Count'])
        for contact, count in contact_counter.items():
            writer.writerow([contact[0], contact[1], count])

    # Counter to aggregate counts for 'Residue 1'
    residue_counter = Counter()

    # Aggregate counts for each unique residue in 'Residue 1' from the summary file
    with open(summary_csv_file, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            residue_counter[row['Residue 1']] += int(row['Count'])

    # Write the additional summary .csv file
    additional_summary_csv_file = os.path.join(output_directory, "residue1_summary.csv")
    with open(additional_summary_csv_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Residue 1', 'Total Count'])
        for residue, total_count in residue_counter.items():
            writer.writerow([residue, total_count])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files to find contacts between specified chains.")
    parser.add_argument("pdb_directory", help="Path to the directory containing .pdb files")
    parser.add_argument("--chain1", default="A", help="Identifier for the first chain (default: A)")
    parser.add_argument("--chain2", default="B", help="Identifier for the second chain (default: B)")
    args = parser.parse_args()

    main(args.pdb_directory, args.chain1, args.chain2)