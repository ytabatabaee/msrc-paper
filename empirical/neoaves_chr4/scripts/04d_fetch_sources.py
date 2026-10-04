#!/usr/bin/env python3
"""Resume pinned chr4avian and historical ASTRAL sources; verify before reuse."""
import json
import subprocess
from datetime import datetime, timezone
from stage4d_core import *

CHR_COMMIT='17b408a622c9589fe9f6ca8a434c150974e39cd4'
ASTRAL_COMMIT='887669674518f38c13c108d006f936cd54d77467'

def main():
    dest=BASE/'external/chr4avian'
    index=json.loads((dest/'repository-index.json').read_text())
    lookup={r['path']:r for r in index['tree']}
    wanted={'README.md':'README.md', 'genetreesupport/README.md':'genetreesupport-README.md',
            'genetreesupport/delta_quartet.R':'delta_quartet.R',
            'genetreesupport/draw-movingaverage.r':'draw-movingaverage.r',
            'genetreesupport/all.stat.xz':'all.stat.xz',
            'genetreesupport/removed-count.tsv':'removed-count.tsv',
            'genetreesupport/taxonsampling.pdf':'taxonsampling.pdf',
            'RYcoding-Lie/README.md':'RYcoding-Lie-README.md',
            'RYcoding-Lie/outlier-indices.txt':'outlier-indices.txt',
            'root-to-tip/README.md':'root-to-tip-README.md'}
    prior={r['filename']:r for r in read_tsv(dest/'download_manifest.tsv')} if (dest/'download_manifest.tsv').exists() else {}
    rows=[]
    for source,name in wanted.items():
        info=lookup[source]; p=dest/name
        url=f'https://raw.githubusercontent.com/smirarab/chr4avian/{CHR_COMMIT}/{source}'
        if name in prior:
            assert sha256(p)==prior[name]['sha256']; rows.append(prior[name]); continue
        if not p.exists() or p.stat().st_size != info['size']:
            subprocess.run(['curl','-fL','-C','-','--retry','5',url,'-o',str(p)],check=True)
        assert p.stat().st_size==info['size']
        # Git blob identity is independent of our local SHA256 inventory.
        import hashlib
        blob=hashlib.sha1(b'blob '+str(p.stat().st_size).encode()+b'\0'+p.read_bytes()).hexdigest()
        assert blob==info['sha'], (source,'Git blob mismatch')
        integrity='size_and_git_blob_verified'
        if name.endswith('.xz'):
            subprocess.run(['xz','-t',str(p)],check=True);integrity+=';xz_pass'
        rows.append(dict(filename=name,source_url=url,archive_path=source,
                         download_timestamp=datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat(),
                         byte_size=p.stat().st_size,sha256=sha256(p),integrity=integrity))
        p.chmod(0o444)
        write_tsv(dest/'download_manifest.tsv',rows)
    write_tsv(dest/'download_manifest.tsv',rows)
    zen=dest/'chr4avian-v1.0.0.zip'
    if zen.exists():
        subprocess.run(['unzip','-t',str(zen)],check=True,stdout=subprocess.DEVNULL)
        rows.append(dict(filename=zen.name,source_url='https://zenodo.org/api/records/10699424/files/smirarab/chr4avian-v1.0.0.zip/content',
                         archive_path='Zenodo record 10699424',download_timestamp=datetime.fromtimestamp(zen.stat().st_mtime,timezone.utc).isoformat(),
                         byte_size=zen.stat().st_size,sha256=sha256(zen),integrity='zip_pass'))
        write_tsv(dest/'download_manifest.tsv',rows);zen.chmod(0o444)
    source=BASE/'external/astral5151/source.tar.gz'
    subprocess.run(['gzip','-t',str(source)],check=True)
    subprocess.run(['tar','-tzf',str(source)],check=True,stdout=subprocess.DEVNULL)
    software=[dict(filename=source.name,source_url=f'https://codeload.github.com/smirarab/ASTRAL/tar.gz/{ASTRAL_COMMIT}',
                   archive_path=ASTRAL_COMMIT,download_timestamp=datetime.fromtimestamp(source.stat().st_mtime,timezone.utc).isoformat(),
                   byte_size=source.stat().st_size,sha256=sha256(source),integrity='gzip_pass;tar_pass')]
    write_tsv(source.parent/'download_manifest.tsv',software);source.chmod(0o444)
    combined=read_tsv(EXTERNAL/'download_manifest.tsv')+rows+software
    write_tsv(RESULTS/'stage4d_download_manifest.tsv',combined)
    (EXTERNAL/'checksums.sha256').write_text(''.join(f"{r['sha256']}  {r['filename']}\n" for r in read_tsv(EXTERNAL/'download_manifest.tsv')))

if __name__=='__main__':main()
