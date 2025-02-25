import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import threading

def split_genome(input_file, output_folder, split_every):
    """Function to split genome into specified chunk sizes"""
    try:
        if not input_file or not output_folder or not split_every:
            log_output("❌ Error: Please provide all required inputs.")
            return

        split_every = int(split_every)

        if not os.path.exists(input_file):
            log_output(f"❌ Error: File not found -> {input_file}")
            return

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Read the input genome file
        with open(input_file, "r") as infile:
            header = infile.readline().strip()  # Read the header line
            sequence = infile.read().replace("\n", "")  # Read full sequence

        # Generate output file name with split size
        base_name = os.path.basename(input_file)
        base_name, ext = os.path.splitext(base_name)
        output_file = os.path.join(output_folder, f"{base_name}_split_{split_every}{ext}")

        with open(output_file, "w") as outfile:
            for i in range(0, len(sequence), split_every):
                start = i + 1
                end = min(i + split_every, len(sequence))
                outfile.write(f"{header}_coord_{start}_{end}\n")
                outfile.write(sequence[i:end] + "\n")

        log_output(f"🎉 Successfully split genome into chunks of {split_every} bases.")
        log_output(f"✅ Output saved to: {output_file}")

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
        title="Select Genome File Path",
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
    split_every = split_entry.get().strip()

    if not input_file or not output_folder or not split_every:
        messagebox.showerror("Error", "Please provide all required inputs.")
        return

    try:
        split_every = int(split_every)
        if split_every <= 0:
            raise ValueError("Split size must be greater than zero.")
    except ValueError:
        messagebox.showerror("Error", "Invalid split size.")
        return

    log_output("🔄 Splitting genome... Please wait.")
    threading.Thread(target=split_genome, args=(input_file, output_folder, split_every), daemon=True).start()

# ----------------------- GUI DESIGN -----------------------
root = tk.Tk()
root.title("Genome Splitter")
root.geometry("700x500")  # Larger window
root.minsize(700, 500)  # Allow resizing

# Frame for better UI organization
frame = tk.Frame(root, padx=20, pady=20)
frame.pack(fill="both", expand=True)

# Input File
tk.Label(frame, text="📂 Genome File Path:", font=("Arial", 12, "bold")).grid(row=0, column=0, sticky="w")
input_entry = tk.Entry(frame, width=60)
input_entry.grid(row=0, column=1, padx=10, pady=5)
tk.Button(frame, text="Browse", command=browse_file).grid(row=0, column=2, padx=5, pady=5)

# Output Folder
tk.Label(frame, text="📁 Output Folder:", font=("Arial", 12, "bold")).grid(row=1, column=0, sticky="w")
output_entry = tk.Entry(frame, width=60)
output_entry.grid(row=1, column=1, padx=10, pady=5)
tk.Button(frame, text="Browse", command=browse_folder).grid(row=1, column=2, padx=5, pady=5)

# Split Size (instead of number of splits)
tk.Label(frame, text="🔢 Bases per Split:", font=("Arial", 12, "bold")).grid(row=2, column=0, sticky="w")
split_entry = tk.Entry(frame, width=10)
split_entry.grid(row=2, column=1, padx=10, pady=5, sticky="w")
split_entry.insert(0, "1000")  # Default value

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
