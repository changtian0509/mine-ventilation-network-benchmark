import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
def patch(path, old, new):
    s = open(path, encoding='utf-8').read()
    n = s.count(old)
    if n == 0:
        print('MISS', path, old[:60]); return
    s = s.replace(old, new)
    open(path, 'w', encoding='utf-8', newline='\n').write(s)
    print('OK', path, n)
patch('refresh_scan.py', '"q_precise":1e-4,"maxCount":maxCount,"h_precise":0.1', '"q_precise":0.01,"maxCount":maxCount,"h_precise":0.1')
patch('refresh_scan_metrics.py', '"q_precise":1e-4,"maxCount":maxCount,"h_precise":0.1', '"q_precise":0.01,"maxCount":maxCount,"h_precise":0.1')
patch('refresh_merit_data.py', 'dg["maxCount"]=800; dg["q_precise"]=1e-2; dg["h_precise"]=5.0', 'dg["maxCount"]=800; dg["q_precise"]=1e-2; dg["h_precise"]=0.1')
