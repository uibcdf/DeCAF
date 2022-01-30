import yaml

def heal(arg):

    output = []

    for ii in arg:
        if type(ii) is list:
            for jj in ii:
                if type(jj) is list:
                    for kk in jj:
                        if type(kk) is list:
                            for ll in kk:
                                if ll not in output:
                                    output.append(ll)
                        else:
                            if kk not in output:
                                output.append(kk)
                else:
                    if jj not in output:
                        output.append(jj)
        else:
            if ii not in output:
                output.append(ii)

    return output

with open('requirements.yaml') as fff:
    all_requirements = yaml.load(fff, Loader=yaml.FullLoader)

## Production

env_dict={}
env_dict["channels"]=heal(all_requirements["production"]["channels"])
env_dict["dependencies"]=heal(all_requirements["production"]["dependencies"])
fff = open("conda-envs/production_env.yaml", "w")
yaml.dump(env_dict, fff, sort_keys=False)
fff.close()

## Development

env_dict={}
env_dict["channels"]=heal(all_requirements["development"]["channels"])
env_dict["dependencies"]=heal(all_requirements["development"]["dependencies"])
fff = open("conda-envs/development_env.yaml", "w")
yaml.dump(env_dict, fff, sort_keys=False)
fff.close()

## Test

env_dict={}
env_dict["channels"]=heal(all_requirements["test"]["channels"])
env_dict["dependencies"]=heal(all_requirements["test"]["dependencies"])
fff = open("conda-envs/test_env.yaml", "w")
yaml.dump(env_dict, fff, sort_keys=False)
fff.close()

## Setup

env_dict={}
env_dict["channels"]=heal(all_requirements["setup"]["channels"])
env_dict["dependencies"]=heal(all_requirements["setup"]["dependencies"])
fff = open("conda-envs/setup_env.yaml", "w")
yaml.dump(env_dict, fff, sort_keys=False)
fff.close()

