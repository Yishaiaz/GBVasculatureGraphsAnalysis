import os
import json
from graph_tool.all import *

"""
exit codes:
1 - unknown error
111 - could not find configuration file 
1112 - could not parse configuration file as json
112 - could not import clear map or compile it
113 - could not load entire graph data 
"""

try:
    with open('job_conf_file.txt', 'r') as f:
        print("FOUND FILE!, the path to the ssd dir is:")
        configuration_file_dict = json.load(f)
    ssd_path_to_dir = configuration_file_dict['ssd_path_to_dir']
except FileNotFoundError as e:
    print("could not find file, but here is the cwd listdir")
    print(f"CWD:{os.getcwd()}")
    print(os.listdir(os.getcwd()))
    exit(111)
except Exception as e:
    print(f'could not parse configuration file json: {e}')
    exit(1112)
print(f"SSD dir path: {ssd_path_to_dir}")

try:
    import sys
    sys.path.insert(0, '/home/yishaiaz/ClearMap2')
    import ClearMap.Compile
except ImportError as e:
    print(f"could not import ClearMap")
    exit(112)

# load graph:

entire_graph_path = os.path.join(ssd_path_to_dir, 'data_graph.gt')

try:
    entire_graph = load_graph(entire_graph_path)
except FileNotFoundError as e:
    print(f"File not found at: {entire_graph_path}")
    exit(113)
except Exception as e:
    print(f"Uknown expection occured: {e}")
    exit(1)

# read vertex annotations
try:
    import ClearMap.Alignment.Annotation as ano
except ImportError as e:
    print(f"could not import annotation module from clearmap {e}")
    exit(1)

try:
    label = entire_graph.vertex_annotation()
    for lvl in range(1, 7):
        label_leveled_converted = ano.convert_label(label, key='order', value='order', level=lvl)
        print(f"level={lvl}\nlabels converted:{label_leveled_converted}\n")
    # vertex_filter = label_leveled == 6
    # gs = entire_graph.sub_graph(vertex_filter=vertex_filter)
except Exception as e:
    print(e)
    exit(1)






