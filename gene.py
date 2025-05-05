import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import seaborn as sns
import numpy as np
from collections import Counter

class ModernGenomicExtractor:
    def __init__(self, root):
        self.root = root
        self.root.title("Modern Genomic Extractor")
        self.root.geometry("800x600")
        
        # Style configuration
        self.style = ttk.Style()
        self.style.configure("Custom.TButton", padding=10, font=('Helvetica', 10))
        self.style.configure("Title.TLabel", font=('Helvetica', 14, 'bold'))
        
        self.file_path = tk.StringVar()
        self.create_gui()
        
    def create_gui(self):
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Title
        title_label = ttk.Label(
            main_frame, 
            text="Genomic Data Extractor", 
            style="Title.TLabel"
        )
        title_label.pack(pady=10)
        
        # File selection frame
        file_frame = ttk.LabelFrame(main_frame, text="File Selection", padding=10)
        file_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Entry(file_frame, textvariable=self.file_path, width=60).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            file_frame, 
            text="Browse", 
            command=self.browse_file, 
            style="Custom.TButton"
        ).pack(side=tk.LEFT, padx=5)
        
        # Action buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(
            button_frame, 
            text="Extract Data", 
            command=self.extract_data, 
            style="Custom.TButton"
        ).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(
            button_frame, 
            text="Save Results", 
            command=self.save_results, 
            style="Custom.TButton"
        ).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(
            button_frame, 
            text="Clear", 
            command=self.clear_results, 
            style="Custom.TButton"
        ).pack(side=tk.LEFT, padx=5)
        
        # Notebook for different views
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Create tabs
        self.source_tab = ttk.Frame(self.notebook)
        self.features_tab = ttk.Frame(self.notebook)
        self.origin_tab = ttk.Frame(self.notebook)
        self.stats_tab = ttk.Frame(self.notebook)
        self.viz_tab = ttk.Frame(self.notebook)
        
        self.notebook.add(self.source_tab, text="Source Info")
        self.notebook.add(self.features_tab, text="Features")
        self.notebook.add(self.origin_tab, text="Origin")
        self.notebook.add(self.stats_tab, text="Statistics")
        self.notebook.add(self.viz_tab, text="Visualizations")
        
        # Create text widgets for each tab
        self.source_text = tk.Text(self.source_tab, height=20, width=70)
        self.features_text = tk.Text(self.features_tab, height=20, width=70)
        self.origin_text = tk.Text(self.origin_tab, height=20, width=70)
        self.stats_text = tk.Text(self.stats_tab, height=20, width=70)
        
        # Add scrollbars
        for text_widget in [self.source_text, self.features_text, 
                          self.origin_text, self.stats_text]:
            scrollbar = ttk.Scrollbar(text_widget.master, orient="vertical", 
                                    command=text_widget.yview)
            text_widget.configure(yscrollcommand=scrollbar.set)
            text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Create visualization buttons
        viz_buttons_frame = ttk.Frame(self.viz_tab)
        viz_buttons_frame.pack(side=tk.TOP, pady=5)
        
        ttk.Button(
            viz_buttons_frame,
            text="Feature Distribution",
            command=self.plot_features,
            style="Custom.TButton"
        ).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(
            viz_buttons_frame,
            text="GC Distribution",
            command=self.plot_gc_distribution,
            style="Custom.TButton"
        ).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(
            viz_buttons_frame,
            text="Base Composition",
            command=self.plot_base_composition,
            style="Custom.TButton"
        ).pack(side=tk.LEFT, padx=5)
        
        # Create frame for plots
        self.plot_frame = ttk.Frame(self.viz_tab)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(fill=tk.X, pady=5)

    def browse_file(self):
        file = filedialog.askopenfilename(
            filetypes=[
                ("GenBank files", "*.gb"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )
        self.file_path.set(file)
        self.status_var.set(f"Selected file: {file}")

    def extract_data(self):
        try:
            file_path = self.file_path.get()
            if not file_path:
                messagebox.showerror("Error", "Please select a file first!")
                return
                
            self.status_var.set("Processing file...")
            record = SeqIO.read(file_path, "genbank")
            self.current_record = record
            
            # Clear all text widgets
            self.clear_results()
            
            # Source Information
            source_info = "SOURCE INFORMATION\n" + "="*50 + "\n"
            source_info += f"Organism: {record.annotations.get('organism', 'N/A')}\n"
            source_info += f"Taxonomy: {', '.join(record.annotations.get('taxonomy', ['N/A']))}\n"
            source_info += f"Accession: {record.id}\n"
            source_info += f"Description: {record.description}\n"
            self.source_text.insert(tk.END, source_info)
            
            # Features
            for feature in record.features:
                feature_info = f"Type: {feature.type}\n"
                feature_info += f"Location: {feature.location}\n"
                feature_info += "Qualifiers:\n"
                for key, value in feature.qualifiers.items():
                    feature_info += f"  {key}: {value}\n"
                feature_info += "-"*50 + "\n"
                self.features_text.insert(tk.END, feature_info)
            
            # Origin
            origin_info = "SEQUENCE\n" + "="*50 + "\n"
            sequence = str(record.seq)
            for i in range(0, len(sequence), 60):
                origin_info += f"{i+1:>8} {sequence[i:i+60]}\n"
            self.origin_text.insert(tk.END, origin_info)
            
            # Statistics
            stats = self.calculate_statistics(record)
            stats_info = "SEQUENCE STATISTICS\n" + "="*50 + "\n"
            stats_info += f"Total Length: {stats['Total Length']} bp\n"
            stats_info += f"GC Content: {stats['GC Content']}%\n"
            stats_info += f"Total Features: {stats['Feature Count']}\n"
            stats_info += f"Molecular Type: {stats['Molecular Type']}\n"
            stats_info += f"Topology: {stats['Topology']}\n\n"
            
            stats_info += "Feature Distribution:\n" + "-"*30 + "\n"
            for feature_type, count in stats['Feature Types'].items():
                stats_info += f"{feature_type}: {count}\n"
            
            self.stats_text.insert(tk.END, stats_info)
            
            self.status_var.set("Data extracted successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to process file: {str(e)}")
            self.status_var.set("Error processing file")

    def calculate_statistics(self, record):
        stats = {
            "Total Length": len(record.seq),
            "GC Content": round(
                (record.seq.count('G') + record.seq.count('C')) / len(record.seq) * 100, 2
            ),
            "Feature Count": len(record.features),
            "Feature Types": {},
            "Molecular Type": record.annotations.get("molecule_type", "Unknown"),
            "Topology": record.annotations.get("topology", "Unknown")
        }
        
        for feature in record.features:
            stats["Feature Types"][feature.type] = stats["Feature Types"].get(feature.type, 0) + 1
            
        return stats

    def plot_features(self):
        if not hasattr(self, 'current_record'):
            messagebox.showerror("Error", "Please extract data first!")
            return
            
        self.clear_plot_frame()
        
        feature_counts = Counter(feature.type for feature in self.current_record.features)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.bar(feature_counts.keys(), feature_counts.values())
        
        ax.set_title('Feature Distribution')
        ax.set_xlabel('Feature Type')
        ax.set_ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}',
                   ha='center', va='bottom')
        
        plt.tight_layout()
        self.show_plot(fig)

    def plot_gc_distribution(self):
        if not hasattr(self, 'current_record'):
            messagebox.showerror("Error", "Please extract data first!")
            return
            
        self.clear_plot_frame()
        
        window_size = 100
        sequence = str(self.current_record.seq)
        gc_content = []
        
        for i in range(0, len(sequence) - window_size, window_size):
            window = sequence[i:i+window_size]
            gc = (window.count('G') + window.count('C')) / window_size * 100
            gc_content.append(gc)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(range(len(gc_content)), gc_content)
        ax.set_title('GC Content Distribution')
        ax.set_xlabel('Position (window size: 100bp)')
        ax.set_ylabel('GC Content (%)')
        
        plt.tight_layout()
        self.show_plot(fig)

    def plot_base_composition(self):
        if not hasattr(self, 'current_record'):
            messagebox.showerror("Error", "Please extract data first!")
            return
            
        self.clear_plot_frame()
        
        sequence = str(self.current_record.seq)
        base_counts = Counter(sequence)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        labels = base_counts.keys()
        sizes = base_counts.values()
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
        ax1.set_title('Base Composition')
        
        bars = ax2.bar(base_counts.keys(), base_counts.values())
        ax2.set_title('Base Counts')
        ax2.set_ylabel('Count')
        
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom')
        
        plt.tight_layout()
        self.show_plot(fig)

    def clear_plot_frame(self):
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

    def show_plot(self, fig):
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def save_results(self):
        save_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if save_path:
            try:
                with open(save_path, 'w') as f:
                    f.write("SOURCE INFORMATION\n\n")
                    f.write(self.source_text.get(1.0, tk.END))
                    f.write("\nFEATURES\n\n")
                    f.write(self.features_text.get(1.0, tk.END))
                    f.write("\nORIGIN\n\n")
                    f.write(self.origin_text.get(1.0, tk.END))
                    f.write("\nSTATISTICS\n\n")
                    f.write(self.stats_text.get(1.0, tk.END))
                self.status_var.set("Results saved successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save results: {str(e)}")

    def clear_results(self):
        self.source_text.delete(1.0, tk.END)
        self.features_text.delete(1.0, tk.END)
        self.origin_text.delete(1.0, tk.END)
        self.stats_text.delete(1.0, tk.END)
        self.status_var.set("Ready")

def main():
    root = tk.Tk()
    app = ModernGenomicExtractor(root)
    root.mainloop()

if __name__ == "__main__":
    main()