# SIMPLER-MAGIC
SIMPLER MAGIC: Synthesis and In-memory MaPping of Logic Execution in a single Row for Memristor Aided loGIC

## Dependencies
In order to use SIMPLER-MAGIC, you will need a Linux machine with:
1. Python 3.6
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
Max_num_gates=20000    ;The maximum number of allowed gates
ROW_SIZE=[5,8,10]
generate_json=True
print_mapping=True
print_warnings=True

```
Change the parameters according to your needs.

2. Run:
```sh
python3 simpler_main.py
```
