import sys,os


if __name__ == '__main__':
    result = sys.popen('samtools depth ~/project_formal/170602_XK/output/XK_result/XK-2T-2/XK-2T-2.dedup.bam | cut -f 3')


    if len(sys.argv) >= 2:
        try:
            result = sum([int(_i) for _i in sys.argv[1:]])
        except ValueError:
            print 'Input unregular args.'


if __name__ == '__main__':
    result = os.system('samtools depth %s | head -n 5| cut -f 3 > temp.file' % sys.argv[1])

    print sum()
