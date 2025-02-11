import tkinter as tk
from tkinter import scrolledtext, messagebox, filedialog
import threading
from contextlib import redirect_stdout

def run_program(params):
    # Replace with your actual function call.
    from asari.main import run_asari
    run_asari(params)

def select_directory(params):
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    directory = filedialog.askdirectory(title="Select a Directory")
    root.destroy()
    params['input'] = directory
    return params

class TextRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, s):
        self.text_widget.insert(tk.END, s)
        self.text_widget.see(tk.END)

    def flush(self):
        pass

def run_program_thread(params, text_widget, continue_button):
    with redirect_stdout(TextRedirector(text_widget)):
        run_program(params)
    continue_button.config(state=tk.NORMAL)

def start_program_gui(params):
    root = tk.Tk()
    root.title("Program Output")

    text_area = scrolledtext.ScrolledText(root, wrap=tk.WORD, height=20, width=80)
    text_area.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

    continue_button = tk.Button(root, text="Continue", state=tk.NORMAL, command=root.destroy)
    continue_button.pack(pady=10)

    threading.Thread(target=run_program_thread, args=(params, text_area, continue_button), daemon=True).start()
    root.mainloop()
    return params

def show_disclaimer():
    accepted = [False]
    def on_accept():
        accepted[0] = True
        disclaimer.destroy()
    disclaimer = tk.Tk()
    disclaimer.title("Disclaimer")
    tk.Label(disclaimer, text="Please read and accept the disclaimer before proceeding:\n\n Asari GUI is Experimental, use at your own risk", wraplength=400, justify="left").pack(padx=20, pady=20)
    tk.Button(disclaimer, text="I Accept", command=on_accept).pack(pady=10)
    disclaimer.mainloop()
    return accepted[0]

def create_ui(data):
    result = {}
    def run_callback():
        nonlocal result
        output = {}
        for key, widget in widgets.items():
            if types[key] is bool:
                output[key] = widget.get()
            else:
                val = widget.get()
                if types[key] is int:
                    try:
                        output[key] = int(val)
                    except ValueError:
                        messagebox.showerror("Error", f"Invalid integer for {key}")
                        return
                elif types[key] is float:
                    try:
                        output[key] = float(val)
                    except ValueError:
                        messagebox.showerror("Error", f"Invalid float for {key}")
                        return
                else:
                    output[key] = val
        result = output
        root.destroy()

    root = tk.Tk()
    root.title("Edit Parameters")

    canvas = tk.Canvas(root)
    vsb = tk.Scrollbar(root, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    frame = tk.Frame(canvas)
    canvas.create_window((0, 0), window=frame, anchor="nw")
    frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

    widgets = {}
    types = {}
    row = 0
    for key, value in data.items():
        tk.Label(frame, text=key).grid(row=row, column=0, padx=5, pady=2, sticky="w")
        if isinstance(value, bool):
            var = tk.BooleanVar(value=value)
            chk = tk.Checkbutton(frame, variable=var)
            chk.grid(row=row, column=1, padx=5, pady=2, sticky="w")
            widgets[key] = var
            types[key] = bool
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            entry = tk.Entry(frame)
            entry.insert(0, str(value))
            entry.grid(row=row, column=1, padx=5, pady=2, sticky="w")
            widgets[key] = entry
            types[key] = type(value)
        elif isinstance(value, (str, type(None))):
            entry = tk.Entry(frame)
            if value is None:
                entry.insert(0, "NONE")
            else:
                entry.insert(0, value)
            entry.grid(row=row, column=1, padx=5, pady=2, sticky="w")
            widgets[key] = entry
            types[key] = str
        row += 1

    run_button = tk.Button(root, text="Continue", command=run_callback)
    run_button.pack(side="bottom", pady=10)
    root.mainloop()
    for k, v in result.items():
        if v == "NONE":
            result[k] = None
    return result

def create_selection_ui(options, data, key):
    result = {}

    def set_value(value):
        nonlocal result
        data[key] = value
        result = data
        root.destroy()

    root = tk.Tk()
    root.title("Select an Option")

    for option in options:
        btn = tk.Button(root, text=option, command=lambda opt=option: set_value(opt))
        btn.pack(pady=2)

    root.mainloop()
    return result

if __name__ == "__main__":
    from asari.default_parameters import PARAMETERS
    from asari.main import SUBCOMMANDS, update_peak_detection_params
    show_disclaimer()
    params = select_directory(PARAMETERS)
    params = create_ui(params)
    params = update_peak_detection_params(params)
    params = create_selection_ui(SUBCOMMANDS, params, "run_gui")
    params = start_program_gui(params)
    
