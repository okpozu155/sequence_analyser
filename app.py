from fastapi import FastAPI, Form, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.responses import StreamingResponse
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from io import BytesIO

app = FastAPI()

# Enable CORS for frontend interaction
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/translate")
async def translate_sequence(sequence: str = Form(...)):
    dna_seq = Seq(sequence.upper())
    if len(dna_seq) % 3 != 0:
        # Option 1: Trim the sequence
        trimmed_seq = dna_seq[:len(dna_seq) - len(dna_seq) % 3]
        protein_seq = trimmed_seq.translate()
        return {
            "warning": "Input sequence length was not a multiple of three and was trimmed.",
            "Trimmed DNA": str(trimmed_seq),
            "Protein": str(protein_seq),
        }

    return {"Protein": str(dna_seq.translate())}


@app.post("/gc_content")
async def gc_content(sequence: str = Form(...)):
    gc_count = sequence.upper().count("G") + sequence.upper().count("C")
    gc_percentage = (gc_count / len(sequence)) * 100
    return {"GC_Content": f"{gc_percentage:.2f}%"}

@app.post("/codon_usage")
async def codon_usage(sequence: str = Form(...)):
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-len(sequence)%3, 3)]
    codon_counts = {codon: codons.count(codon) for codon in set(codons)}
    return codon_counts

@app.post("/visualize_gc")
async def visualize_gc(sequence: str = Form(...)):
    # Calculate GC content
    gc_content = [(sequence[:i].count("G") + sequence[:i].count("C")) / i for i in range(1, len(sequence) + 1)]

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(sequence) + 1), gc_content, label="GC Content")
    plt.xlabel("Position in Sequence")
    plt.ylabel("GC Content")
    plt.title("GC Content Across Sequence")
    plt.legend()
    plt.grid()

    # Save the plot to a BytesIO buffer
    buffer = BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)  # Rewind the buffer to the beginning
    plt.close()

    # Return the image as a StreamingResponse
    return StreamingResponse(buffer, media_type="image/png")