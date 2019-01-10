# SIMPLE-MAGIC
SIMPLER MAGIC: Synthesis and In-memory MaPping of Logic Execution in a single Row for Memristor Aided loGIC

## Dependencies
In order to use SIMPLER-MAGIC, you will need a Linux machine with:
1. Python 2.7
2. [ABC Synthesis Tool](https://bitbucket.org/alanmi/abc)

## Manual
1. Configure: in the file simple_conf.cfg you will find the following content:
```ini
[input_output]
input_path=full_adder_1bit.v
; input_format can get one of the values: verilog, blif
input_format=verilog
output_path=full_adder_1bit_output

[abc]
abc_dir_path=/home/adi/abc/alanmi-abc-eac02745facf

[SIMPLER_Mapping]
MAX_D=20000
ROW_SIZE=[5,8,10]

```
Change the parameters according to your needs.

2. Run:
```sh
python simpler_main.py
```