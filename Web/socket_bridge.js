// Socket Bridge - 连接网页和C++程序的中间层
// 使用Node.js运行此脚本: node socket_bridge.js

const http = require('http');
const net = require('net');

const HTTP_PORT = 8765;  // 网页连接的端口
const CPP_PORT = 8765;   // C++程序监听的端口
const CPP_HOST = 'localhost';

console.log('启动Socket桥接服务器...');

// 创建HTTP服务器供网页连接
const httpServer = http.createServer((req, res) => {
    // 处理CORS
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'POST, OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') {
        res.writeHead(200);
        res.end();
        return;
    }

    if (req.method === 'POST' && req.url === '/convert') {
        let body = '';
        
        req.on('data', chunk => {
            body += chunk.toString();
        });

        req.on('end', () => {
            console.log('\n收到网页请求，数据长度:', body.length);
            
            // 连接到C++程序
            const client = new net.Socket();
            let cppResponse = '';

            client.connect(CPP_PORT, CPP_HOST, () => {
                console.log('已连接到C++程序');
                client.write(body);
            });

            client.on('data', (data) => {
                cppResponse += data.toString();
            });

            client.on('end', () => {
                console.log('收到C++响应，数据长度:', cppResponse.length);
                res.writeHead(200, { 'Content-Type': 'application/json' });
                res.end(cppResponse);
                client.destroy();
            });

            client.on('error', (err) => {
                console.error('连接C++程序失败:', err.message);
                res.writeHead(500, { 'Content-Type': 'application/json' });
                res.end(JSON.stringify({ 
                    status: 'error', 
                    message: '无法连接到C++程序，请确保程序正在运行' 
                }));
                client.destroy();
            });
        });
    } else {
        res.writeHead(404);
        res.end('Not Found');
    }
});

httpServer.listen(HTTP_PORT, () => {
    console.log(`HTTP服务器已启动，监听端口: ${HTTP_PORT}`);
    console.log(`等待网页连接...`);
    console.log(`连接到C++程序端口: ${CPP_PORT}`);
    console.log('\n提示：');
    console.log('1. 确保C++程序正在运行');
    console.log('2. 在浏览器中打开 Molecular3D.html');
    console.log('3. 点击"3D转换"按钮进行转换\n');
});
