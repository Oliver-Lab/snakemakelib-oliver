def relative_symlink(input, output):
    """ Helper function to easily symlink two files """
    import os
    from snakemake import shell
    odir = os.path.dirname(output)
    oname = os.path.basename(output)
    relative_path = './' + os.path.relpath(str(input), odir)
    shell("cd {odir}; ln -s {rp} {output}".format(odir=odir, rp=relative_path, output=oname))


def file_md5sum(input):
    import hashlib
    with open(input, 'rb') as IN:
        md5 = hashlib.md5(IN.read()).hexdigest()
    return md5
