# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for kmer manipulation.
'''

# DEPENDENCIES =================================================================

from Bio.SeqIO.FastaIO import SimpleFastaParser
from ggc.args import check_threads
import ggc.fasta as fasta
from heapq import merge
from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import os
import re
import resource
import shutil
import tempfile
from tqdm import tqdm

# CONSTANTS ====================================================================

tqdm.monitor_interval = 0
handler_limit = resource.getrlimit(resource.RLIMIT_NOFILE)[0]

# FUNCTIONS ====================================================================

def rc(seq, ab):
    '''
    Args:
        seq (string): nucleic acid sequence.
        ab (list): alphabet list, with letters and their reverse.

    Return:
        string: reverse complement of seq.
    '''

    assert type([]) == type(ab)
    assert 2 == len(ab)
    assert all([type("") == type(x) for x in ab])

    rab = ab[1].strip().lower()
    ab = ab[0].strip().lower()

    # Check provided string
    seq = seq.lower()
    for c in seq:
        if not c in ab:
            print('ERROR: provided string conflicts with the selected alphabet.')
            return

    # Calculate reverse
    r = seq[::-1]

    # Calculate reverse complement
    rc = r
    for ci in range(len(ab)):
        rc.replace(ab[ci], rab[ci])

    return(rc.toupper())

def parallel_sets_union(plist, opath, threads = 1, progress = True,
    frec = None):
    '''Make union of sorted sets of kmers written to disk, in parallel.

    Args:
        plist (list): list of tuples with abs.paths to sorted kmer sets
                      and their size.
        opath (str): path to output kmer set after union.
        threads (int): number of threads for parallelization.
        progress (bool): show progress bar.
        frec (fun): intermediate generator function to run on FASTA records.

    Returns:
        int: number of output sequences
        str: output path.
    '''

    threads = check_threads(threads) # Check thread count
    verbose = 11 if progress else 0     # Parallel verbosity

    # Calculate sorting batch size
    batch_size = min(int(len(plist) / float(threads)), handler_limit - 24)

    # Parallelized union
    if batch_size >= len(plist):
        plist = [sets_union(plist, opath, progress)]
    else:
        print(" Performing parallel sorted unique kmer set union...")
        print(" Batch size: %d" % batch_size)
        plist = Parallel(n_jobs = threads, verbose = verbose)(
        delayed(sets_union)(plist[i:min(len(plist), i + batch_size)],
            progress = False, headless = True, frec = frec)
            for i in range(0, len(plist), batch_size))

    # If multiple outputs from parallelization, perform union
    if 1 != len(plist):
        out = sets_union(plist, opath, progress, frec = frec)
        for (p, n) in plist:
            if p != opath: os.remove(p)
        return(out)
    else:
        plist = plist[0]
        if plist[0] != opath:
            shutil.copyfile(plist[0], opath)
            os.remove(plist[0])
        return((opath, plist[1]))

def sets_union(plist, opath = None, progress = False, headless = False,
    frec = None):
    '''Make union of sorted sets of kmers written to disk.

    First merge kmer sets with a k-way merge (heapq.merge) to maintain the
    sorting, then iterate through the file and write only the first occurrence
    of each kmer in the output.

    Args:
        plist (list): list of tuples with abs.paths to sorted kmer sets
                      and their size.
        opath (str): path to output kmer set after union.
        progress (bool): show progress bar.
        frec (fun): intermediate generator function to run on FASTA records.

    Returns:
        int: number of output sequences
        str: output path.
    '''

    def iterate_sorted_kset(path):
        '''Iterates over sequences in header-less FASTA file.

        Args:
            path (str): path to header-less FASTA.
        '''

        # Kill if file not found
        if not os.path.isfile(path): return

        # Prepare iterator
        with open(path, 'r') as IH:
            for record in SimpleFastaParser(IH):
                yield(record)

                    
    def write_seq(OH, record, c):
        '''Writes one sequence to a FASTA file.

        Args:
            OH (TextIOWrapper): handle to output file.
            record (tuple): (header, sequence).
            c (int): counter for id.
        '''
        rid = "candidate_%d" % c
        rdescr = " ".join(record[0].split(" ")[1:])
        OH.write(">%s %s\n%s\n" % (rid, rdescr, record[1]))
        return(c + 1)

    def do_pass(): return(curr_c == 1 and 0 != len(curr_rec[1]))

    # Log status
    if progress: print(" Performing sorted unique kmer set union...")

    if type(None) == type(opath):
        (H, opath) = tempfile.mkstemp(".setUnion_output.fa"); os.close(H)

    # Perform union in batches to avoid handle soft limit
    (H, tpath) = tempfile.mkstemp(".setUnion_fasta.fa"); os.close(H)
    for i in range(0, len(plist), handler_limit - 24):
        tplist = [plist[j][0] for j in range(i,
            min(len(plist), i + handler_limit - 24))]
        tplist.append(tpath)

        # Prepare list of iterable generators
        generators = [iterate_sorted_kset(p) for p in tplist]

        # Merge iterables
        crawl = merge(*generators, key = lambda x: x[1])
        crawl = frec(crawl) if type(lambda x: x) == type(frec) else crawl

        # Write as output those sequences that appear only once
        curr_rec = ["", ""]   # Current sequence
        curr_c = 1      # Occurences counter
        out_c = 1       # Candidate counter
        with open(opath, 'w+') as OH:
            nseq = sum([t[1] for t in plist])
            if progress: crawl = tqdm(crawl, total = nseq)
            
            for record in crawl:
                # Increase occurrence counter
                if record[1] == curr_rec[1]: curr_c += 1
                else:
                    # Write out previous seq
                    if do_pass(): out_c = write_seq(OH, curr_rec, out_c)

                    # Move to next
                    curr_rec = record
                    curr_c = 1

            # Write out last seq
            if do_pass(): out_c = write_seq(OH, curr_rec, out_c)

        if not plist[-1][0] in tplist:
            os.rename(opath, tpath)
    os.remove(tpath)

    return((opath, out_c - 1))

def uniq_as_fasta(k, record, klim, threads = 1, allow_non_ACTUG = False,
    encountered = None, single = False, fasta_delim = "=", frec = None):
    '''Extract unique kmers (k-characters substring) from a sequence (string).
    Converts the input sequence in a FASTA-like file with records of up to klim
    kmers, and uses the uniq_fasta function to avoid keeping the whole kmer set
    in memory.

    Args:
        k (int): oliqonucleotide (substring) length in nt.
        record (tuple): (header, sequence).
        klim (int): upper limit of kmers in a batch, must be >= 100.
        threads (int): number of threads for parallelization.
        allow_non_ACTUG (bool): whether to allow non-ACTUG characters.
        encountered (set): set of already encountered kmers to skip.
        single (bool): single record run.
        fasta_delim (str): FASTA header key-val delimiter.
        frec (fun): intermediate generator function to run on FASTA records.

    Returns:
        list: sorted unique kmer list.
        str: path to output file.
    '''

    threads = check_threads(threads) # Check thread count
    if klim < 100: klim = 100  # Set minimum batch size

    seq = record[1]
    chrom = re.findall(r'pos%s(.+):(.+)-(.+);' % fasta_delim, record[0])
    cstart = 0 if 0 == len(chrom) else int(chrom[0][1])
    cend = 0 if 0 == len(chrom) else int(chrom[0][2])
    chrom = record[0].split(" ")[0] if 0 == len(chrom) else chrom[0][0]

    # Convert to temporary FASTA file
    (H, tipath) = tempfile.mkstemp(".uaf.fa"); os.close(H)
    c = 0
    with open(tipath, "w+") as TIH:
        for tstart in range(0, len(seq), klim):
            c += 1
            tend = min(len(seq), tstart + klim + k - 1)
            rid = record[0].split(" ")[0]

            # Build and write FASTA record
            rid = "seq_%d" % c
            rdescr = "pos%s%s:%d-%d;" % (fasta_delim, chrom,
                tstart + 1 + cstart - 1, tend + 1 + cend - 1)
            TIH.write(">%s %s\n%s\n" % (rid, rdescr, seq[tstart:tend]))

    # Set verbosity
    progress = False if 1 == threads else True

    # Generate unique kmers
    (H, topath) = tempfile.mkstemp(".uaf_kset.fa"); os.close(H)
    nseq = uniq_fasta(k, tipath, topath, threads = threads,
        allow_non_ACTUG = allow_non_ACTUG, encountered = encountered,
        progress = progress, fasta_delim = fasta_delim)
    os.remove(tipath)

    if not single:
        print(" Extracting kmers from batch output...")
        with open(topath, "r") as IH:
            # Prepare sequence generator
            def fastagen(IH):
                for r in SimpleFastaParser(IH): yield(r)

            nrecs = fasta.count_records(topath)
            if type(lambda x: x) == type(frec):
                crawler = frec(fastagen(IH))
            else:
                crawler = fastagen(IH)
            if progress: crawler = tqdm(crawler, total = nrecs)

            kset = []
            for r in crawler: kset.append((r[1], r[0]))
        os.remove(topath)
        return((kset, topath))
    else:
        return(([], topath))

def uniq_fasta(k, fpath, opath, klim = 0, threads = 1, allow_non_ACTUG = False,
    encountered = None, progress = False, fasta_delim = "=", frec = None,
    reverseAB = None, **kwargs):
    '''Extract unique kmers from a FASTA file.
    Stores the list of uniqued kmers in memory for one FASTA record at a time.

    For each of the M FASTA records, extract unique kmers and write the sorted
    kmer set to disk in a tmp folder. Then, perform M-way merge and iterated
    through the merged sorted file and retain only kmers that appear once.

    Args:
        k (int): oliqonucleotide (substring) length in nt.
        fpath (str): path to input FASTA file.
        opath (str): path to output FASTA file.
        klim (int): upper limit of kmers in a batch, triggers batched analysis.
        threads (int): number of threads for parallelization.
        allow_non_ACTUG (bool): whether to allow non-ACTUG characters.
        encountered (set): set of already encountered kmers to skip.
        progress (bool): show progress bar.
        fasta_delim (str): FASTA header key-val delimiter.
        frec (fun): intermediate generator function to run on FASTA records.

    Returns:
        int: size of sorted unique kmer set written to disk.
    '''

    threads = check_threads(threads) # Check thread count

    # Make temporary folder
    tmpdir = tempfile.mkdtemp("_uniq_fasta")
    tplist = []

    # Open input fasta buffer
    IH = open(fpath, "r")

    # Check if only one record is present in the FASTA file
    single = True if 1 == fasta.count_records(fpath) else False

    # Common arguments
    pargs = {'k':k, 'outdir':tmpdir, 'single':single,
        'fasta_delim':fasta_delim, 'frec':frec}

    def recGenBase(IH):
        for record in SimpleFastaParser(IH):
            yield(record)
    if not type(None) == type(reverseAB):
        def recGen(IH):
            for record in recGenBase(IH):
                yield(record)
                record[0] = f"rc_{record[0]}"
                record[1] = rc(record[1], reverseAB)
                yield(record)
    else:
        recGen = recGenBase


    # Parse FASTA records
    if 1 == threads or 0 != klim:
        # Either normally or in batches (and parallel)
        for record in recGen(IH):
            tplist.append(uniq_record(record,
                klim = klim, threads = threads, progress = progress, **pargs))
    else:
        verbose = 11 if progress else 0 # Parallel verbosity

        # In parallel
        tplist = Parallel(n_jobs = threads, verbose = verbose)(
            delayed(uniq_record)(record, progress = False, **pargs)
            for record in recGen(IH))

    # Remove None objects
    tplist = [p for p in tplist if type(None) != type(p)]

    # Close input buffer
    IH.close()

    if 1 != len(tplist):
        # Perform sorted sets union and make output FASTA
        opath, nseq = parallel_sets_union(tplist, opath,
            threads = threads, progress = progress)
    else:
        # Convert single header-less FASTA to FASTA
        tpath, nseq = tplist[0]

        print(" Moving output...")
        shutil.copyfile(tpath, opath)
        os.remove(tpath)

    # Delete temporary folder
    shutil.rmtree(tmpdir)

    return(nseq)

def uniq_record(record, k, outdir, klim = 0, threads = 1,
    allow_non_ACTUG = False, encountered = None, single = False,
    progress = False, fasta_delim = "=", frec = None):
    '''Extract unique kmers from a single FASTA record.
    Either convert to FASTA for batched analysis, or run normally.

    Args:
        record (tuple): (header, sequence).
        k (int): oliqonucleotide (substring) length in nt.
        outdir (str): path to output folder (usually a temporary one).
        klim (int): upper limit of kmers in a batch, triggers batched analysis.
        threads (int): number of threads for parallelization.
        allow_non_ACTUG (bool): whether to allow non-ACTUG characters.
        encountered (set): set of already encountered kmers to skip.
        single (bool): single record run.
        progress (bool): show progress bar.
        fasta_delim (str): FASTA header key-val delimiter.
        frec (fun): intermediate generator function to run on FASTA records.

    Returns:
        str: path to output file.
        int: number of sequences in output file.
        None: if no kmers are accepted.
    '''

    threads = check_threads(threads) # Check thread count

    # Log record
    if progress: print(" Parsing kmers from record '%s'..." % record[0])

    # Common arguments
    pargs = {'k':k, 'fasta_delim':fasta_delim, 'record':record,
        'allow_non_ACTUG':allow_non_ACTUG, 'encountered':encountered}

    # Extract unique kmers
    if klim > 0:
        # Batched analysis
        kset, topath = uniq_as_fasta(klim = klim,
            threads = threads, single = single, **pargs)

        if single: return((topath, fasta.count_records(topath)))
    else:
        kset, discarded = uniq_seq(progress = progress, **pargs)

        # Sort
        if progress: print(
            " Sorting unique kmer set for record '%s'..." % record[0])
        kset = sorted(kset.items(), key = lambda x: x[0])

    if 0 != len(kset):
        # Write FASTA output
        if progress:
            print(" Writing sorted kmer set for record '%s'..." % record[0])
            
        rid = record[0].split(" ")[0]
        (H, tpath) = tempfile.mkstemp(".txt", "%s." % rid, outdir)
        os.close(H)

        with open(tpath, 'w+') as OH:
            def crawler(kset):
                for (seq, head) in kset: yield((head, seq))

            crawl = crawler(kset)
            crawl = tqdm(crawl, total = len(kset)) if progress else crawl
            crawl = frec(crawl) if type(lambda x: x) == type(frec) else crawl
            for (head, seq) in crawl:
                OH.write(">%s\n%s\n" % (head, seq))
        
        return((tpath, len(kset)))

def uniq_seq(k, record, allow_non_ACTUG = False, encountered = None,
    progress = False, fasta_delim = "=", frec = None):
    '''Extract unique kmers (k-characters substring) from a sequence (string).
    Needs to be able store the list of uniqued kmers in memory.

    Args:
        k (int): oliqonucleotide (substring) length in nt.
        record (tuple): (header, sequence).
        allow_non_ACTUG (bool): whether to allow non-ACTUG characters.
        encountered (set): set of already encountered kmers to skip.
        progress (bool): show progress bar.
        fasta_delim (str): FASTA header key-val delimiter.
        frec (fun): intermediate generator function to run on FASTA records.
    '''

    def update_sets(kmer, name, chrom, start, end):
        '''Update current unique kmer and discarded kmer sets.'''
        
        # Add if never encountered, remove otherwise
        if kmer not in encountered:
            if kmer not in kset.keys():
                header = ">%s pos%s%s:%d-%d" % (
                    name, fasta_delim, chrom, start, end)
                if type(lambda x: x) == type(frec):
                    (header, kmer) = frec((header, kmer))
                kset[kmer] = header[1:]
            else:
                del kset[kmer]
                encountered.add(kmer)

    def kpgen(crawler, seq):
        '''Generate (kmer, position) tuples.'''
        for i in crawler:
            yield((seq[i:(i + k)].upper(), i))

    seq = record[1]
    chrom = re.findall(r'pos%s(.+):(.+)-(.+);' % fasta_delim, record[0])
    cstart = 0 if 0 == len(chrom) else int(chrom[0][1])
    cend = 0 if 0 == len(chrom) else int(chrom[0][2])
    chrom = record[0].split(" ")[0] if 0 == len(chrom) else chrom[0][0]

    # Prepare sets
    kset = {}
    if type(None) == type(encountered):
        encountered = set()

    # Setup crawler
    crawler = range(len(seq) - k + 1)
    crawler = tqdm(crawler) if progress else crawler
    crawler = kpgen(crawler, seq)

    # Extract kmers
    for (kmer, i) in crawler:
        # Update sets if non-ACTUG characters are allowed
        if allow_non_ACTUG: update_sets(kmer, record[0].split(" ")[0],
            chrom, i + cstart - 1, i + k + cstart - 1)
        else:
            # Check for non-ACTUG
            skip = False
            for c in "RYKMSWBDHVN":
                if c in kmer:
                    skip = True
                    break

            # Skip or update
            if not skip: update_sets(kmer, record[0].split(" ")[0],
                chrom, i + cstart - 1, i + k + cstart - 1)

    return((kset, encountered))

# END ==========================================================================

################################################################################
