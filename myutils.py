import re
from typing import Iterable, Tuple


def grab(file: str, start: str, end: str):
    header, body, footer = "", "", ""
    ishead = True
    isfoot = False
    with open(file) as f:
        line = f.readline()
        while line:
            if re.match(start, line):
                ishead = False
            if re.match(end, line):
                isfoot = True

            if ishead:
                header += line
            elif isfoot:
                footer += line
            else:
                body += line

            line = f.readline()

    return (header, body, footer)


def continuous(sequence: Iterable[int], skip: int=1) -> Tuple[Tuple[int], Tuple[int]]:

    lseq = list(sequence)
    seq = iter(sequence)
    continuous = list()
    second = next(seq)
    bpoint = 0
    while True:
        try:
            first = second
            second = next(seq)
            if first + skip != second:
                continuous.append(tuple(lseq[bpoint:lseq.index(second)]))
                bpoint = lseq.index(second)
        except StopIteration:
            continuous.append(tuple(lseq[bpoint:]))
            break
    return tuple(continuous)


if __name__ == "__main__":
    s1 = list(range(2, 16))
    s2 = list(range(30, 42))
    s1.append(25)
    s1.extend(s2)

    print(continuous(s1))
