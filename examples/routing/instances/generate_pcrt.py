#make pcrt benchmarks in the style of Alexander Nadel's FMCAD 16 RUC paper
import random
import argparse

def make_pcrt(f, N, M, C, seed):
    print("Creating benchmark: N %d, M %d, C %d"%(N,M,C))


    width = N*M
    height = N*M
    assert(width>1 or height>1)
    if seed is not None:
        print("Setting seed to %d"%(seed))
        random.seed(seed)
    f.write("G %d %d\n"%(width,height))
    used = set()

    def newPoint(used,shareXY=None):
        if shareXY is None:
            x1 = random.randint(0, width - 1)
            y1 = random.randint(0, height - 1)
            while ((x1, y1) in used):
                x1 = random.randint(0, width - 1)
                y1 = random.randint(0, height - 1)
        else:
            #Generate a new point that shares either x or y with shareXY
            assert(len(shareXY)==2)
            shareX = random.choice([True, False])
            if shareX:
                x1 = shareXY[0]
                y1 = random.randint(0, height - 1)
            else:
                x1 = random.randint(0, width - 1)
                y1 = shareXY[1]
            while ((x1, y1) in used):
                if shareX:
                    x1 = shareXY[0]
                    y1 = random.randint(0, height - 1)
                else:
                    x1 = random.randint(0, width - 1)
                    y1 = shareXY[1]

        used.add((x1, y1))
        return (x1,y1)

    def toInt(x,y):
        assert(x>=0)
        assert(y>=0)
        assert (x < width)
        assert (y < height)
        return y*width + x

    for i in range(N):

        x1,y1 = newPoint(used)
        #make sure the two points are not identical
        x2,y2 = newPoint(used)
        f.write("N %d %d\n"%(toInt(x1,y1),toInt(x2,y2)))

    V = (N*M)^2
    for i in range((C*V)//100):
        x1,y1 = newPoint(used)
        x2, y2 = newPoint(used,(x1,y1))
        f.write("C %d %d\n" % (toInt(x1, y1), toInt(x2, y2)))

    f.close()
    print("Done")

if __name__ == '__main__':
    import sys

    parser = argparse.ArgumentParser(description='Generate PCRT benchmarks, as described in the Alexander Nadel FMCAD16 RUC paper')

    parser.add_argument('N', metavar='N', default = 20, type=int,
                        help='N')
    parser.add_argument('M', metavar='M', type=int,
                        help='M')
    parser.add_argument('constraints', metavar='c', type=int,
                        help='constraints')
    parser.add_argument('seed', type=int, default=None,
                        help='Random seed (default: randomize)')
    parser.add_argument('filename', type=argparse.FileType('w'))

    args = parser.parse_args()

    make_pcrt(args.filename,args.N, args.M, args.constraints, args.seed)
