import yaml
import numpy as np
import os


def _pprint_dict(d, indent=0):
    """Pretty print a dictionary

    Print a dictionary where nesting is visualized with indentation.
    Useful for validation and debugging.

    Args:
        d (dict): The dictionary to print
        indent (int, optional): The amount of indentation to use for this level, default=0

    """
    for key, val in d.items():
        for k in range(indent):
            print("|   ", end="")

        if isinstance(val, dict):
            print(f"{key}:")
            indent += 1
            _pprint_dict(val, indent)
            indent -= 1

        else:
            print(f"{key}: {val}")#" ({type(val)})")


def build_simple_dict(d, parameters_filename, first_call=True, nested=0, path=[]):
    
    if first_call:
        with open(parameters_filename, "w") as fp:
            fp.write("A New Parameters File\n")
            fp.write("=====================\n")
    
    for key, val in d.items():
        spacer = ""
        if nested > 1:
            for k in range(nested-1):
                if k == 0:
                    spacer += "\|\-\-"
                elif k < nested-1:
                    spacer += " \|\-\-"
                
            spacer += "> "

        path.append(key)
        label = ":".join(path[1:]) # 1: don't print parent object, e.g., "name" not "general:name"
        label = f"{spacer}``{label}``"

        if "properties" in val:
            if nested == 0:
                with open(parameters_filename, "a") as fp:
                    fp.write(f"\n{key}\n")
                    fp.write(f"------------------------------------\n\n")

                    parent_description = val.get("description")
                    if parent_description is not None:
                        fp.write(f"{parent_description}\n\n")

                    fp.write(f".. list-table::\n")
                    fp.write(f"  :widths: 20 15 15 10 40\n")
                    fp.write(f"  :header-rows: 1\n\n")

                    fp.write(f"  * - Name\n")
                    fp.write(f"    - Type\n")
                    fp.write(f"    - Units\n")
                    fp.write(f"    - Default\n")
                    fp.write(f"    - Description\n")
                    
            else:
                with open(parameters_filename, "a") as fp:
#                     fp.write(f"  * - \n")
#                     fp.write(f"    - \n")
#                     fp.write(f"    - \n")
#                     fp.write(f"    - \n")
#                     fp.write(f"    - \n")

                    fp.write(f"  * - {label}\n")
                    fp.write(f"    - \n")
                    fp.write(f"    - \n")
                    fp.write(f"    - \n")
                    fp.write(f"    - \n")
                    
            nested, path = build_simple_dict(val["properties"], 
                                       parameters_filename, 
                                       first_call=False, 
                                       nested=nested+1, 
                                       path=path)
            nested -= 1
            
#             if nested == 1:
#                 with open(parameters_filename, "a") as fp:

#                     fp.write(f"  * - \n")
#                     fp.write(f"    - \n")
#                     fp.write(f"    - \n")
#                     fp.write(f"    - \n")
#                     fp.write(f"    - \n")

            
        else:
            with open(parameters_filename, "a") as fp:
                fp.write(f"  * - {label}\n")
                fp.write(f"    - {val['type']}\n")
                fp.write(f"    - NA\n")
                fp.write(f"    - {val['default']}\n")
                fp.write(f"    - {val['description']}\n")
                
        path.pop(-1)
                
    return nested, path


def build_params_rst_from_schema(path_to_schema_file, path_to_params_rst_file):
    with open(path_to_schema_file, "r") as fp:
        schema_dict = yaml.safe_load(fp)

    _, _ = build_simple_dict(schema_dict["properties"], path_to_params_rst_file)


if __name__ == "__main__":

    # Note that these paths should be relative to the root WindSE
    # directory, since that's where the .readthedocs.yaml file is located
    path_to_schema_file = "windse/input_schema.yaml"
    path_to_params_rst_file="doc/source/parameters.rst"

    build_params_rst_from_schema(path_to_schema_file, path_to_params_rst_file)
