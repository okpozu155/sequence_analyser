from fastapi import FastAPI, Form, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.responses import StreamingResponse
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
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
    try:
        # Clean the input sequence: remove invalid characters
        clean_sequence = ''.join([char for char in sequence.upper() if char in "ATGC"])
        
        # Check if the cleaned sequence is empty
        if not clean_sequence:
            return {"error": "No valid DNA sequence provided after cleaning input."}
        
        # Trim the sequence to a multiple of 3
        trimmed_sequence = clean_sequence[:len(clean_sequence) - len(clean_sequence) % 3]
        
        # Translate the sequence
        dna_seq = Seq(trimmed_sequence)
        protein_seq = dna_seq.translate()
        
        return {
            "Original DNA": sequence,
            "Cleaned DNA": clean_sequence,
            "Trimmed DNA": trimmed_sequence,
            "Protein": str(protein_seq),
        }
    except TranslationError as e:
        # Handle translation errors gracefully
        return {"error": f"Translation failed: {str(e)}"}
    except Exception as e:
        # Catch unexpected errors
        return {"error": f"An unexpected error occurred: {str(e)}"}


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