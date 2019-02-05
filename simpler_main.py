import ConfigParser
import os
import tempfile
import SIMPLER_Mapping
import ast

def main():

    # Read configuration parameters
    config = ConfigParser.ConfigParser()
    config.readfp(open('simpler_conf.cfg'))
    input_path = config.get('input_output', 'input_path')
    input_format = config.get('input_output', 'input_format')
    abc_dir_path = config.get('abc', 'abc_dir_path')
    #BenchmarkStrings = ast.literal_eval(config.get("SIMPLER_Mapping", "BenchmarkStrings"))
    Max_num_gates = config.get('SIMPLER_Mapping', 'Max_num_gates')
    ROW_SIZE = [int(i) for i in ast.literal_eval(config.get("SIMPLER_Mapping", "ROW_SIZE"))]
    output_path = config.get('input_output', 'output_path')
    generate_json = config.get('SIMPLER_Mapping', 'generate_json')
    print_mapping = config.get('SIMPLER_Mapping', 'print_mapping')
    print_warnings = config.get('SIMPLER_Mapping', 'print_warnings')

    abc_exe_path = os.path.join(abc_dir_path, "abc")
    abc_rc_path = os.path.join(abc_dir_path, "abc.rc")

    # Create abc script
    abc_script = file('abc_script_template.abc', 'rb').read()
    abc_script = abc_script.replace('abc_rc_path', abc_rc_path)
    abc_script = abc_script.replace('input.blif', input_path)
    if input_format == 'verilog':
        abc_script = abc_script.replace('read_blif', 'read_verilog')
    abc_script = abc_script.replace('lib.genlib', 'mcnc1_nor2.genlib')
    abc_output_path = tempfile.mktemp()
    abc_script = abc_script.replace('output.v', abc_output_path)

    # Run abc script
    abc_script_path = tempfile.mktemp()
    file(abc_script_path, "wb").write(abc_script)
    os.system('%s -f "%s"' % (abc_exe_path, abc_script_path))

    # Mapping into the memory array
    SIMPLER_Mapping.SIMPLER_Main([abc_output_path], Max_num_gates, ROW_SIZE, input_path.split(".")[0], generate_json, print_mapping, print_warnings)
    

    # Clean files
    os.remove(abc_script_path)
    # os.remove(abc_output_path)

if __name__ == "__main__":
    main()
