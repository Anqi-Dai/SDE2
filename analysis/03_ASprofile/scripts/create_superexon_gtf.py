import re
fn = '/usr4/bs831/adai/bubhub-home/SDE2/reference/v27/gencode.v27.annotation.gtf'

fields = (
        'gene_id','transcript_id','gene_name','gene_type',
        'transcript_type'
    )

def print_synthetic_transcript_lines(curr_transcript) :
    # saw a new transcript, print out the super transcript and
    # exon records
    synth_transcript_id = '"{}_super"'.format(curr_transcript['gene_id'].replace('"',''))
    curr_transcript['transcript_id'] = synth_transcript_id
    curr_transcript['transcript_type'] = '"synthetic_super_exon"'
    attr_str = []
    for f in fields :
        attr_str.append('{} {}'.format(f,curr_transcript[f]))

    # print out the synthetic transcript record
    print('\t'.join([
        curr_transcript['chrm'],
        curr_transcript['src'],
        'transcript',
        str(curr_transcript['start']),
        str(curr_transcript['end']),
        '.',
        curr_transcript['strand'],
        '.',
        '; '.join(attr_str)
    ]))

    # print out super exon
    print('\t'.join([
        curr_transcript['chrm'],
        curr_transcript['src'],
        'exon',
        str(curr_transcript['start']),
        str(curr_transcript['end']),
        '.',
        curr_transcript['strand'],
        '.',
        '; '.join(attr_str)
    ]))

comments = []
with open(fn) as f :
    curr_transcript = {}
    for line in f :
        if not line.startswith('#') :
            parts = line.split('\t')
            chrm, src, feature, start, end, phase, strand, score, attrs = parts

            attrs_d = {}
            for attr in attrs.split(';') :
                attr = attr.strip()
                if len(attr) != 0 :
                    key, val = attr.split(' ',1)
                    attrs_d[key] = val

            if curr_transcript.get('gene_id') is None :
                curr_transcript = attrs_d
                curr_transcript['chrm'] = chrm
                curr_transcript['src'] = src
                curr_transcript['start'] = int(start)
                curr_transcript['end'] = int(end)
                curr_transcript['strand'] = strand
            elif curr_transcript['gene_id'] == attrs_d['gene_id'] :
                curr_transcript['start'] = min(curr_transcript['start'],int(start))
                curr_transcript['end'] = max(curr_transcript['end'],int(end))
            else :

                print_synthetic_transcript_lines(curr_transcript)

                # recording a new curr_transcript
                curr_transcript = attrs_d
                curr_transcript['chrm'] = chrm
                curr_transcript['src'] = src
                curr_transcript['start'] = int(start)
                curr_transcript['end'] = int(end)
                curr_transcript['strand'] = strand

        print(line.strip())

    # this handles the last gene in the file
    curr_transcript = attrs_d
    curr_transcript['chrm'] = chrm
    curr_transcript['src'] = src
    curr_transcript['start'] = int(start)
    curr_transcript['end'] = int(end)
    curr_transcript['strand'] = strand

    print_synthetic_transcript_lines(curr_transcript)
