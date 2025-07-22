from importlib.resources import files

def get_data_file_path(filename):
    return str(files("cathapult.data").joinpath(filename))