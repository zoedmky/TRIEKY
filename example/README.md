# Example Datasets

This folder contains toy VCF records and sample configuration files to test TRIEKY.

## DMPK Test

Run the following command from this directory:

```bash
python ../trieky.py targets_example.tsv dmpk_example.vcf dmpk_output.txt
```
(Or using uv: uv run python ../trieky.py targets_example.tsv dmpk_example.vcf dmpk_output.txt)

The generated dmpk_output.txt report should be identical to dmpk_expected.txt.

## RFC1 Test
Run the following command from this directory:

```bash
python ../trieky.py targets_example.tsv rfc1_example.vcf rfc1_output.txt
```
(Or using uv: uv run python ../trieky.py targets_example.tsv rfc1_example.vcf rfc1_output.txt)

The generated rfc1_output.txt report should be identical to rfc1_expected.txt.
