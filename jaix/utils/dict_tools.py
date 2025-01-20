# https://stackoverflow.com/questions/14692690/access-nested-dictionary-items-via-a-list-of-keys


def nested_get(dic, keys):
    for key in keys:
        dic = dic[key]
    return dic


def nested_set(dic, keys, value):
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value


def nested_del(dic, keys):
    for key in keys[:-1]:
        dic = dic[key]
    del dic[keys[-1]]
