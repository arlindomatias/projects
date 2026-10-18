from transformers import pipeline, AutoTokenizer, AutoModelForSeq2SeqLM

text = (
    "Transformers are a powerful deep learning architecture designed "
    "for processing sequential data, such as natural language, by "
    "utilizing a \"self-attention\" mechanism to evaluate the relevance "
    "of all other tokens in a sequence simultaneously."
)

# --- Summarization (BART) ---
sum_tokenizer = AutoTokenizer.from_pretrained("facebook/bart-large-cnn")
sum_model = AutoModelForSeq2SeqLM.from_pretrained("facebook/bart-large-cnn")

inputs = sum_tokenizer(text, return_tensors="pt", truncation=True, max_length=1024)
summary_ids = sum_model.generate(**inputs, max_length=30, min_length=5, do_sample=False)
summary = sum_tokenizer.decode(summary_ids[0], skip_special_tokens=True)

print(summary)

# --- Sentiment classification ---
classifier = pipeline("sentiment-analysis", top_k=None)  # or "text-classification"

classification = classifier(text)

print(classification)

# --- Translation (NLLB) ---
tr_model_name = "facebook/nllb-200-distilled-600M"
tr_tokenizer = AutoTokenizer.from_pretrained(tr_model_name, src_lang="eng_Latn")
tr_model = AutoModelForSeq2SeqLM.from_pretrained(tr_model_name)

tr_inputs = tr_tokenizer(text, return_tensors="pt")
tr_ids = tr_model.generate(
    **tr_inputs,
    forced_bos_token_id=tr_tokenizer.convert_tokens_to_ids("por_Latn"),
    max_new_tokens=200,
)
translation = tr_tokenizer.batch_decode(tr_ids, skip_special_tokens=True)[0]

print(translation)