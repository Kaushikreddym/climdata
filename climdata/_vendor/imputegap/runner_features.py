from imputegap.recovery.explainer import Explainer
from imputegap.recovery.manager import TimeSeries
from imputegap.tools import utils

dataset = "temperature"

ts = TimeSeries()
ts.load_series(data=utils.search_path(dataset), header=False)

categories, features, _ = Explainer.load_configuration()
characteristics, descriptions = Explainer.extractor_pycatch(data=ts.data, features_categories=categories, features_list=features, do_catch24=False)

p = "./dataset/docs/"+dataset+"/features_"+dataset+".txt"
with open(p, 'w') as f:
    for desc in descriptions:
        key, category, description = desc
        if key in characteristics:
            value = characteristics[key]
            f.write(f"|{category}|{description}|{value}|\n")
        else:
            f.write(f"Warning: Key '{key}' not found in characteristics!\n")
print(f"Table exported to {p}")
