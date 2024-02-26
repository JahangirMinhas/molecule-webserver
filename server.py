from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import parse_qs
import _molecule

class MyServer(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == '/':
            self.send_response(200)
            self.send_header('Content-type', 'text/html')
            self.end_headers()
            self.wfile.write('''<html>
 <head>
 <title> File Upload </title>
 </head>
 <body>
 <h1> File Upload </h1>
 <form action="molecule" enctype="multipart/form-data" method="post">
 <p>
 <input type="file" id="sdf_file" name="filename"/>
 </p>
 <p>
 <input type="submit" value="Upload"/>
 </p>
 </form>
 </body>
</html>'''.encode())
        else:
            self.send_error(404, 'File Not Found')

    def do_POST(self):
        if self.path == '/molecule':
            content_length = int(self.headers['Content-Length'])
            post_data = self.rfile.read(content_length)
            form_data = parse_qs(post_data)
            if b'filename' in form_data:
                file_content = form_data[b'filename'][0]
                molecule = Molecule.from_sdf(file_content)
                svg_data = molecule.to_svg()
                self.send_response(200)
                self.send_header('Content-type', 'image/svg+xml')
                self.end_headers()
                self.wfile.write(svg_data)
            else:
                self.send_error(400, 'Bad Request')
        else:
            self.send_error(404, 'File Not Found')

def run(server_class=HTTPServer, handler_class=MyServer, port=5000):
    server_address = ('', port)
    httpd = server_class(server_address, handler_class)
    print(f'Starting server on port {port}...')
    httpd.serve_forever()

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        port = int(sys.argv[1])
        run(port=port)
    else:
        print('Usage: python server.py <port>')