from transformers import pipeline

summarizer = pipeline("summarization", model="facebook/bart-large-cnn")

texto = 'Transformers are a powerful deep learning architecture designed for processing sequential data, such as natural language, by utilizing a "self-attention" mechanism to evaluate the relevance of all other tokens in a sequence simultaneously'

resumo = summarizer(texto, max_length=30,
                    min_length=5,
                    do_sample=False)

print(resumo)

classifier = pipeline("sentiment-analysis", top_k=None)  # ou "text-classification"

classification = classifier(texto)

print(classification)

translator = pipeline(
    "translation",
    model="facebook/nllb-200-distilled-600M",
    src_lang="eng_Latn",
    tgt_lang="por_Latn",
)

translation = translator(texto)
print(translation[0]["translation_text"])