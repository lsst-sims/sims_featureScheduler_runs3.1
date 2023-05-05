import numpy as np

if __name__ == "__main__":
    weights = [12,24,50,75,100,125,150]
    with open('ug.sh', 'w') as file:
        for w1 in weights:
            for w2 in weights:
                file.write('python ug_y1.py --utw %i --gtw %i \n' % (w1, w2), )
