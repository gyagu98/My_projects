import os
import time
import threading
from Bio import SeqIO
import tkinter as tk
from tkinter import messagebox
from tkinter.ttk import Progressbar
from Cocoa import NSOpenPanel, NSSavePanel


def split_genome(input_fasta, output_dir, num_chunks=None, chunk_size=None, progress_callback=None):
    """
    Splits a genome into either a fixed number of chunks or fixed-size chunks.

    :param input_fasta: Path to the input FASTA file
    :param output_dir: Directory to store output FASTA files
    :param num_chunks: Number of chunks to split the genome into
    :param chunk_size: Size of each chunk in base pairs (e.g., 20000 for 20 KB)
    :param progress_callback: Function to update progress bar
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    full_sequence = "".join(str(record.seq) for record in records)  # Combine all sequences into one
    total_length = len(full_sequence)
    
    if num_chunks:
        chunk_size = total_length // num_chunks
    elif chunk_size:
        num_chunks = (total_length + chunk_size - 1) // chunk_size  # Ensure all bases are included
    else:
        raise ValueError("Either num_chunks or chunk_size must be provided.")
    
    for i in range(num_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, total_length)
        chunk_seq = full_sequence[start:end]
        
        chunk_filename = os.path.join(output_dir, f"genome_chunk_{i+1}.fasta")
        with open(chunk_filename, "w") as out_file:
            out_file.write(f">chunk_{i+1}\n{chunk_seq}\n")
        
        if progress_callback:
            progress_callback(i + 1, num_chunks)
        
        print(f"Saved: {chunk_filename}")
        time.sleep(0.5)  # Simulate processing time for progress update
    
    print("Splitting complete!")


class GenomeSplitterGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Genome Splitter GUI")
        self.root.geometry("500x300")
        
        tk.Label(root, text="Input FASTA File:").pack()
        self.input_file_entry = tk.Entry(root, width=50)
        self.input_file_entry.pack()
        tk.Button(root, text="Browse", command=self.browse_input_file).pack()
        
        tk.Label(root, text="Output Directory:").pack()
        self.output_dir_entry = tk.Entry(root, width=50)
        self.output_dir_entry.pack()
        tk.Button(root, text="Browse", command=self.browse_output_dir).pack()
        
        tk.Label(root, text="Number of Chunks (leave blank for size-based splitting):").pack()
        self.num_chunks_entry = tk.Entry(root, width=10)
        self.num_chunks_entry.pack()
        
        tk.Label(root, text="Chunk Size in Base Pairs (leave blank if using chunks):").pack()
        self.chunk_size_entry = tk.Entry(root, width=10)
        self.chunk_size_entry.pack()
        
        self.progress_bar = Progressbar(root, orient="horizontal", length=300, mode="determinate")
        self.progress_bar.pack(pady=10)
        
        tk.Button(root, text="Split Genome", command=self.run_split).pack()
        tk.Button(root, text="Exit", command=root.quit).pack()
    
    def browse_input_file(self):
        panel = NSOpenPanel.openPanel()
        panel.setCanChooseFiles_(True)
        panel.setCanChooseDirectories_(False)
        panel.setAllowsMultipleSelection_(False)
        
        if panel.runModal():
            filename = panel.URLs()[0].path()
            self.input_file_entry.delete(0, tk.END)
            self.input_file_entry.insert(0, filename)
    
    def browse_output_dir(self):
        panel = NSSavePanel.savePanel()
        panel.setCanCreateDirectories_(True)
        panel.setCanChooseDirectories_(True)
        panel.setAllowsMultipleSelection_(False)
        
        if panel.runModal():
            dirname = panel.URL().path()
            self.output_dir_entry.delete(0, tk.END)
            self.output_dir_entry.insert(0, dirname)
    
    def update_progress(self, current, total):
        percent = int((current / total) * 100)
        self.progress_bar["value"] = percent
        self.root.update_idletasks()
    
    def run_split(self):
        input_fasta = self.input_file_entry.get().strip()
        output_dir = self.output_dir_entry.get().strip()
        num_chunks = self.num_chunks_entry.get().strip()
        chunk_size = self.chunk_size_entry.get().strip()
        
        num_chunks = int(num_chunks) if num_chunks.isdigit() else None
        chunk_size = int(chunk_size) if chunk_size.isdigit() else None
        
        if not input_fasta or not output_dir:
            messagebox.showerror("Error", "Please select both input file and output directory")
            return
        
        threading.Thread(target=split_genome, args=(input_fasta, output_dir, num_chunks, chunk_size, self.update_progress), daemon=True).start()


if __name__ == "__main__":
    root = tk.Tk()
    app = GenomeSplitterGUI(root)
    root.mainloop()
