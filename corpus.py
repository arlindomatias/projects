import tkinter as tk
from tkinter import filedialog, messagebox
import fitz  # PyMuPDF
import os

def process_pdf():
    # abrir seletor de arquivo
    file_path = filedialog.askopenfilename(
        title="Selecione um PDF",
        filetypes=[("PDF Files", "*.pdf")]
    )
    if not file_path:
        return

    # abrir PDF
    doc = fitz.open(file_path)

    # criar pasta para salvar resultados
    output_dir = os.path.splitext(file_path)[0] + "_output"
    os.makedirs(output_dir, exist_ok=True)

    # pegar metadados
    metadata = doc.metadata
    meta_text = ""
    for k, v in metadata.items():
        if v:
            meta_text += f"{k}: {v}\n"

    # salvar texto
    text_file = os.path.join(output_dir, "texto.txt")
    with open(text_file, "w", encoding="utf-8") as f:
        for page in doc:
            f.write(page.get_text() + "\n")

    # salvar imagens
    image_dir = os.path.join(output_dir, "figuras")
    os.makedirs(image_dir, exist_ok=True)
    img_count = 0
    for page_num, page in enumerate(doc, start=1):
        for img in page.get_images(full=True):
            xref = img[0]
            pix = fitz.Pixmap(doc, xref)
            img_count += 1
            img_path = os.path.join(image_dir, f"figura_{page_num}_{img_count}.png")
            if pix.n < 5:  # RGB ou GRAY
                pix.save(img_path)
            else:  # CMYK, converter para RGB
                pix = fitz.Pixmap(fitz.csRGB, pix)
                pix.save(img_path)
            pix = None

    # exibir mensagem
    messagebox.showinfo(
        "Processo concluído",
        f"Metadados extraídos:\n{meta_text}\n\n"
        f"Texto salvo em: {text_file}\n"
        f"{img_count} imagens salvas em: {image_dir}"
    )

# criar interface
root = tk.Tk()
root.title("Extrator de Papers PDF")

btn = tk.Button(root, text="Selecionar PDF", command=process_pdf)
btn.pack(pady=20, padx=20)

root.mainloop()

