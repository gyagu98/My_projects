import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import threading

def split_genome(input_file, output_folder, num_splits):
    """Function to split genome file into parts"""
    try:
        if not input_file or not output_folder or not num_splits:
            log_output("❌ Error: Please provide all required inputs.")
            return
        
        num_splits = int(num_splits)

        if not os.path.exists(input_file):
            log_output(f"❌ Error: File not found -> {input_file}")
            return

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Read the genome file
        with open(input_file, "r") as f:
            lines = f.readlines()

        # Find sequence lines (ignore headers for now)
        sequences = [line.strip() for line in lines if not line.startswith(">")]

        # Combine sequence into one long string
        full_sequence = "".join(sequences)

        # Split the genome into chunks
        chunk_size = len(full_sequence) // num_splits
        output_files = []

        for i in range(num_splits):
            start = i * chunk_size
            end = (i + 1) * chunk_size if i < num_splits - 1 else len(full_sequence)
            chunk = full_sequence[start:end]
            
            output_filename = os.path.join(output_folder, f"split_{i+1}.fasta")
            with open(output_filename, "w") as out_f:
                out_f.write(f">split_{i+1}\n")
                for j in range(0, len(chunk), 80):  # Formatting to 80-char lines
                    out_f.write(chunk[j:j+80] + "\n")

            output_files.append(output_filename)
            log_output(f"✅ Saved chunk {i+1} to {output_filename}")

        log_output(f"🎉 Successfully split genome into {num_splits} parts!")

    except Exception as e:
        log_output(f"❌ Error: {str(e)}")

def log_output(message):
    """Logs messages in the GUI console"""
    log_text.config(state=tk.NORMAL)
    log_text.insert(tk.END, message + "\n")
    log_text.config(state=tk.DISABLED)
    log_text.yview(tk.END)

def browse_file():
    """Open file dialog to select input genome file"""
    filename = filedialog.askopenfilename(
        title="Select Genome File",
        initialdir=os.getcwd(),
        filetypes=[("FASTA Files", "*.fasta *.fa *.fna"), ("All Files", "*.*")]
    )
    if filename:
        input_entry.delete(0, tk.END)
        input_entry.insert(0, filename)

def browse_folder():
    """Open folder dialog to select output directory"""
    foldername = filedialog.askdirectory(
        title="Select Output Folder",
        initialdir=os.getcwd()
    )
    if foldername:
        output_entry.delete(0, tk.END)
        output_entry.insert(0, foldername)

def start_splitting():
    """Start genome splitting process in a separate thread"""
    input_file = input_entry.get().strip()
    output_folder = output_entry.get().strip()
    num_splits = split_entry.get().strip()

    if not input_file or not output_folder or not num_splits:
        messagebox.showerror("Error", "Please provide all required inputs.")
        return

    try:
        num_splits = int(num_splits)
        if num_splits <= 0:
            raise ValueError("Number of splits must be greater than zero.")
    except ValueError:
        messagebox.showerror("Error", "Invalid number of splits.")
        return

    log_output("🔄 Splitting genome... Please wait.")
    threading.Thread(target=split_genome, args=(input_file, output_folder, num_splits), daemon=True).start()

# ----------------------- GUI DESIGN -----------------------
root = tk.Tk()
root.title("Genome Splitter")
root.geometry("700x500")  # Larger window
root.minsize(700, 500)  # Allow resizing

# Frame for better UI organization
frame = tk.Frame(root, padx=20, pady=20)
frame.pack(fill="both", expand=True)

# Input File
tk.Label(frame, text="📂 Genome File:", font=("Arial", 12, "bold")).grid(row=0, column=0, sticky="w")
input_entry = tk.Entry(frame, width=60)
input_entry.grid(row=0, column=1, padx=10, pady=5)
tk.Button(frame, text="Browse", command=browse_file).grid(row=0, column=2, padx=5, pady=5)

# Output Folder
tk.Label(frame, text="📁 Output Folder:", font=("Arial", 12, "bold")).grid(row=1, column=0, sticky="w")
output_entry = tk.Entry(frame, width=60)
output_entry.grid(row=1, column=1, padx=10, pady=5)
tk.Button(frame, text="Browse", command=browse_folder).grid(row=1, column=2, padx=5, pady=5)

# Number of Splits
tk.Label(frame, text="🔢 Number of Splits:", font=("Arial", 12, "bold")).grid(row=2, column=0, sticky="w")
split_entry = tk.Entry(frame, width=10)
split_entry.grid(row=2, column=1, padx=10, pady=5, sticky="w")
split_entry.insert(0, "3")  # Default value

# Start Button - Larger and Centered
start_btn = tk.Button(
    frame, text="🚀 Start Splitting", font=("Arial", 14, "bold"), bg="#4CAF50", fg="white", command=start_splitting
)
start_btn.grid(row=3, column=0, columnspan=3, pady=15, sticky="ew")

# Log Output Area
tk.Label(frame, text="📜 Log Output:", font=("Arial", 12, "bold")).grid(row=4, column=0, sticky="w")
log_text = tk.Text(frame, height=12, width=80, state=tk.DISABLED, bg="black", fg="white", font=("Courier", 10))
log_text.grid(row=5, column=0, columnspan=3, padx=10, pady=5)

# Allow resizing
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

root.mainloop()
