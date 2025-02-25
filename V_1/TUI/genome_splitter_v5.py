import os
import time
from Bio import SeqIO
from textual.app import App, ComposeResult
from textual.containers import Vertical
from textual.widgets import Button, Input, Label, ProgressBar
from textual.worker import Worker, WorkerState

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

class GenomeSplitterApp(App):
    """A simple TUI app for genome splitting with a working progress bar and exit option."""

    def compose(self) -> ComposeResult:
        yield Label("Input FASTA File:")
        self.input_file = Input()
        yield self.input_file
        
        yield Label("Output Directory:")
        self.output_dir = Input()
        yield self.output_dir
        
        yield Label("Number of Chunks (leave blank for size-based splitting):")
        self.num_chunks = Input()
        yield self.num_chunks
        
        yield Label("Chunk Size in Base Pairs (leave blank if using chunks):")
        self.chunk_size = Input()
        yield self.chunk_size
        
        yield Button("Split Genome", id="split")
        yield Button("Exit", id="exit")
        
        self.progress_bar = ProgressBar(total=100)
        yield self.progress_bar
    
    def update_progress(self, current, total):
        percent = int((current / total) * 100)
        self.progress_bar.progress = percent
        self.refresh()
    
    def run_split(self):
        self.run_worker(self.worker_split, thread=True, exclusive=True)
    
    def worker_split(self):
        input_fasta = self.input_file.value.strip()
        output_dir = self.output_dir.value.strip()
        num_chunks = int(self.num_chunks.value.strip()) if self.num_chunks.value.strip().isdigit() else None
        chunk_size = int(self.chunk_size.value.strip()) if self.chunk_size.value.strip().isdigit() else None
        
        split_genome(input_fasta, output_dir, num_chunks=num_chunks, chunk_size=chunk_size, progress_callback=self.update_progress)
    
    def on_button_pressed(self, event):
        if event.button.id == "split":
            self.run_split()
        elif event.button.id == "exit":
            self.exit()

if __name__ == "__main__":
    app = GenomeSplitterApp()
    app.run()
