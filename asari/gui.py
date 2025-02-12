import tkinter as tk
from tkinter import scrolledtext, messagebox, filedialog
import threading
from contextlib import redirect_stdout, redirect_stderr
from asari.default_parameters import PARAMETERS
from asari.main import SUBCOMMANDS, update_peak_detection_params, run_asari

def run_program(params):
    run_asari(params)

def select_directory(params):
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(title="Select Input Directory")
    params['input'] = directory
    root.destroy()
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
    redirector = TextRedirector(text_widget)
    continue_button.config(state=tk.DISABLED)
    with redirect_stdout(redirector), redirect_stderr(redirector):
        run_program(params)
    continue_button.config(state=tk.NORMAL)


def start_program_gui(params):
    root = tk.Tk()
    root.title("Asari Processing")

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
    def on_decline():
        exit()

    disclaimer = tk.Tk()
    disclaimer.title("Disclaimer")

    DISCLAIMER = """
    
    ** Asari GUI is Experimental **

    Asari is a tool for LC and GC-MS metabolomics data preprocessing. 
    However, the GUI is still in the experimental stage and may not be fully functional. 
    Please use the GUI with caution and report any issues to us on GitHub.

    By clicking "I Accept", you agree to the terms and conditions of the disclaimer.

    """
    tk.Label(disclaimer, text=DISCLAIMER, wraplength=400, justify="center").pack(padx=20, pady=20)
    tk.Button(disclaimer, text="I Accept", command=on_accept).pack(pady=10)
    tk.Button(disclaimer, text="I Decline", command=on_decline).pack(pady=10)

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
        try:
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
        except:
            pass

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

def main_gui():
    if not show_disclaimer():
        return
    params = select_directory(PARAMETERS)
    if params['input'] is None:
        messagebox.showerror("Error", "No directory selected. Exiting.")
        return
    params = create_ui(params)
    if 'autoheight' not in params:
        messagebox.showerror("Error", "Parameters issue (code 1). Exiting.")
        return
    try:
        params = update_peak_detection_params(params)
        params = create_selection_ui(SUBCOMMANDS, params, "run_gui")
    except:
        messagebox.showerror("Error", "Parameters issue (code 2). Exiting.")

    try:
        params = start_program_gui(params)
    except Exception as e:
        messagebox.showerror("Thanks for trying Asari GUI: ", e)

if __name__ == "__main__":
    main_gui()
    
