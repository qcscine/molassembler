from pymol import cmd, stored

def fit_all_to_first():
  all_objects_list = cmd.get_object_list(selection='(all)')
  for name in all_objects_list[1:]:
    cmd.fit(name, all_objects_list[0])

cmd.extend("fit_all_to_first", fit_all_to_first);
