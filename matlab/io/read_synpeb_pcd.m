function pc = read_synpeb_pcd(filename)
% Read SynPEB PCD (float32 scalar fields), including LZF binary_compressed.
% Base MATLAB only. File ordering is preserved: point index=(row-1)*width+col.
fid = fopen(filename,'rb','ieee-le');
assert(fid>=0,'Cannot open %s',filename);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
hdr = struct();
while true
    line = fgetl(fid);
    assert(ischar(line),'Truncated PCD header: %s',filename);
    line = strtrim(line);
    if isempty(line) || line(1)=='#', continue; end
    t = strsplit(line); key = upper(t{1});
    hdr.(key) = t(2:end);
    if strcmp(key,'DATA'), break; end
end
n = str2double(hdr.POINTS{1});
nf = numel(hdr.FIELDS);
assert(all(str2double(hdr.SIZE)==4) && all(strcmp(hdr.TYPE,'F')) ...
    && all(str2double(hdr.COUNT)==1), ...
    'This reader expects SynPEB scalar float32 fields.');
switch hdr.DATA{1}
    case 'binary_compressed'
        sizes = fread(fid,2,'uint32=>double');
        assert(numel(sizes)==2 && sizes(2)==4*n*nf,'Invalid payload size.');
        src = fread(fid,sizes(1),'*uint8');
        assert(numel(src)==sizes(1),'Truncated compressed payload.');
        raw = lzf_decode(src,sizes(2));
        fields = reshape(typecast(raw,'single'),n,nf);
    case 'binary'
        fields = fread(fid,[nf n],'single=>single')';
    case 'ascii'
        fields = fscanf(fid,'%f',[nf n])';
    otherwise
        error('Unsupported PCD encoding: %s',hdr.DATA{1});
end
assert(size(fields,1)==n,'Incorrect PCD point count.');
ix = find(strcmp(hdr.FIELDS,'x')); iy=find(strcmp(hdr.FIELDS,'y'));
iz = find(strcmp(hdr.FIELDS,'z')); il=find(strcmp(hdr.FIELDS,'intensity'));
assert(~isempty(il),'Missing intensity reference labels.');
pc.xyz = double(fields(:,[ix iy iz]));
pc.labels = double(fields(:,il));
pc.width = str2double(hdr.WIDTH{1});
pc.height = str2double(hdr.HEIGHT{1});
pc.filename = filename;
end

function out = lzf_decode(src,nout)
out = zeros(nout,1,'uint8'); ip=1; op=1; ns=numel(src);
while ip<=ns
    ctrl=double(src(ip)); ip=ip+1;
    if ctrl<32
        len=ctrl+1;
        assert(ip+len-1<=ns && op+len-1<=nout,'Invalid LZF literal.');
        out(op:op+len-1)=src(ip:ip+len-1); ip=ip+len; op=op+len;
    else
        len=floor(ctrl/32); ref=op-mod(ctrl,32)*256-1;
        if len==7, len=len+double(src(ip)); ip=ip+1; end
        ref=ref-double(src(ip)); ip=ip+1; len=len+2;
        assert(ref>=1 && op+len-1<=nout,'Invalid LZF back-reference.');
        % LZF permits overlap, so copy sequentially in short chunks.
        remaining=len;
        while remaining>0
            chunk=min(remaining,op-ref);
            out(op:op+chunk-1)=out(ref:ref+chunk-1);
            op=op+chunk; ref=ref+chunk; remaining=remaining-chunk;
        end
    end
end
assert(op==nout+1,'LZF decompressed length mismatch.');
end

