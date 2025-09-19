import tkinter as tk 
from tkinter import filedialog

def select_file():
    root = tk.Tk()
    root.withdraw()  # hide the root window
    file_path = filedialog.askopenfilename(
        title="Select a file",
        filetypes=[("All files", "*.*"), ("Text files", "*.txt"), ("CSV files", "*.csv")]
    )
    return file_path

select_file()