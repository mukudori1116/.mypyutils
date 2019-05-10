import re


def grab(file: str, start: str, end: str):
    with open(file) as f:
        lines = f.readlines()
    ishead = True
    isfoot = False
    header, footer = "", ""
    body = []
    for line in lines:
        if re.match(end, line):
            isfoot = True

        if ishead:
            header += line
        elif isfoot:
            footer += line
        else:
            body.append(line)

        if re.match(start, line):
            ishead = False
    return (header, "".join(body), footer)
