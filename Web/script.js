// script.js
const colors = {
    C: '#FFFFFF',
    H: '#A9A9A9',
    O: '#87CEFA',
    N: '#FFD700',
    S: '#FFB6C1'
};

let atoms = [];
let bonds = [];
let selectedAtoms = new Set();
let selectedBonds = new Set();
let currentTool = 'C';
let currentBond = 'single';
let nextAtomId = 0;
let nextBondId = 0;
let mode = 'draw';
let dragStartPos = { x: 0, y: 0 };
let dragEndPos = { x: 0, y: 0 };
let isDraggingAtom = false;
let dragTargetId = null;
let dragStartOffset = { x: 0, y: 0 };
let isDraggingSelection = false;

// 确保在 HTML 文件中定义了 svg 元素
const svg = document.getElementById('svg');
const canvas = document.getElementById('canvas');

let isSelecting = false;
let startX = 0;
let startY = 0;
let selectionRect;
let selectionBox = { left: 0, top: 0, right: 0, bottom: 0 };


const ATOM_RADIUS = 10;
const ATOM_FONT_SIZE = '14px';


const BOND_STROKE_WIDTH = 2 //键的宽度
const BOND_END_OFFSET = 10; // 键缩短的长度
const BOND_TYPE_OFFSET = 2;
const WEDGE_WIDTH = 8;

const HIGHLIGHT_COLOR = 'rgba(255, 255, 200, 0.5)';
const HOVER_COLOR = 'rgba(255, 255, 200, 0.4)';
const CLICK_THRESHOLD = 10;


function draw() {
    svg.innerHTML = '';
    const atomMap = new Map(atoms.map(a => [a.id, a]));
    const bondCount = atoms.reduce((acc, atom) => {
        acc[atom.id] = bonds.filter(b => b.a1 === atom.id || b.a2 === atom.id).length;
        return acc;
    }, {});

    const createLine = (x1, y1, x2, y2, strokeWidth = BOND_STROKE_WIDTH, strokeColor = colors.C) => {
        const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttribute('x1', x1);
        line.setAttribute('y1', y1);
        line.setAttribute('x2', x2);
        line.setAttribute('y2', y2);
        line.setAttribute('stroke', strokeColor);
        line.setAttribute('stroke-width', strokeWidth);
        svg.appendChild(line);
    };

    const createWedge = (x1, y1, x2, y2, type, fillColor = colors.C) => {
        const dx = x2 - x1;
        const dy = y2 - y1;
        const dist = Math.sqrt(dx * dx + dy * dy);
        const angle = Math.atan2(dy, dx);
        const perpAngle = angle + Math.PI / 2;

        if (type === 'wedge-up') {
            const offset = WEDGE_WIDTH / 2;
            const endOffX = offset * Math.cos(perpAngle);
            const endOffY = offset * Math.sin(perpAngle);

            const p1 = `${x1},${y1}`;
            const p2 = `${x2 + endOffX},${y2 + endOffY}`;
            const p3 = `${x2 - endOffX},${y2 - endOffY}`;

            const polygon = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
            polygon.setAttribute('points', `${p1} ${p2} ${p3}`);
            polygon.setAttribute('fill', fillColor);
            polygon.setAttribute('stroke', 'none');
            svg.appendChild(polygon);
        } else if (type === 'wedge-down') {
            const dashLength = 3;
            const dashGap = 2;
            let currentDist = 0;

            while (currentDist < dist) {
                const startRatio = currentDist / dist;
                const endRatio = Math.min(1, (currentDist + dashLength) / dist);

                const startWidth = WEDGE_WIDTH * startRatio;
                const endWidth = WEDGE_WIDTH * endRatio;

                const perpOffset1_start = { x: (startWidth / 2) * Math.cos(perpAngle), y: (startWidth / 2) * Math.sin(perpAngle) };
                const perpOffset1_end = { x: (endWidth / 2) * Math.cos(perpAngle), y: (endWidth / 2) * Math.sin(perpAngle) };

                const perpOffset2_start = { x: -(startWidth / 2) * Math.cos(perpAngle), y: -(startWidth / 2) * Math.sin(perpAngle) };
                const perpOffset2_end = { x: -(endWidth / 2) * Math.cos(perpAngle), y: -(endWidth / 2) * Math.sin(perpAngle) };

                const seg_x1 = x1 + startRatio * dx + perpOffset1_start.x;
                const seg_y1 = y1 + startRatio * dy + perpOffset1_start.y;
                const seg_x2 = x1 + endRatio * dx + perpOffset1_end.x;
                const seg_y2 = y1 + endRatio * dy + perpOffset1_end.y;

                const seg_x3 = x1 + startRatio * dx + perpOffset2_start.x;
                const seg_y3 = y1 + startRatio * dy + perpOffset2_start.y;
                const seg_x4 = x1 + endRatio * dx + perpOffset2_end.x;
                const seg_y4 = y1 + endRatio * dy + perpOffset2_end.y;

                const polygon = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
                polygon.setAttribute('points', `${seg_x1},${seg_y1} ${seg_x2},${seg_y2} ${seg_x4},${seg_y4} ${seg_x3},${seg_y3}`);
                polygon.setAttribute('fill', fillColor);
                polygon.setAttribute('stroke', 'none');
                svg.appendChild(polygon);

                currentDist += dashLength + dashGap;
            }
        }
    };

    bonds.forEach(b => {
        const a1 = atomMap.get(b.a1);
        const a2 = atomMap.get(b.a2);
        if (!a1 || !a2) return;

        const dx_total = a2.x - a1.x;
        const dy_total = a2.y - a1.y;
        const dist_total = Math.sqrt(dx_total * dx_total + dy_total * dy_total);
        const angle = Math.atan2(dy_total, dx_total);
        const perpAngle = angle + Math.PI / 2;

        let startOffset = 0;
        let endOffset = 0;
        if (a1.elem !== 'C') {
            startOffset = BOND_END_OFFSET;
        }
        if (a2.elem !== 'C') {
            endOffset = BOND_END_OFFSET;
        }
        const dist = dist_total - startOffset - endOffset;
        if (dist <= 0) return;

        const startX = a1.x + startOffset * Math.cos(angle);
        const startY = a1.y + startOffset * Math.sin(angle);
        const endX = a2.x - endOffset * Math.cos(angle);
        const endY = a2.y - endOffset * Math.sin(angle);

        const strokeColor = selectedBonds.has(b.id) ? 'rgb(255, 255, 0)' : colors.C; // 键的高亮颜色

        // 绘制点击热区
        const hotspot = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        hotspot.setAttribute('x1', startX);
        hotspot.setAttribute('y1', startY);
        hotspot.setAttribute('x2', endX);
        hotspot.setAttribute('y2', endY);
        hotspot.setAttribute('class', `selectable-bond${selectedBonds.has(b.id) ? ' selected' : ''}`);
        hotspot.setAttribute('data-bond-id', b.id);
        svg.appendChild(hotspot);

        const offset = BOND_TYPE_OFFSET;

        if (b.type === 'double') {
            const offX = offset * Math.cos(perpAngle);
            const offY = offset * Math.sin(perpAngle);
            createLine(startX + offX, startY + offY, endX + offX, endY + offY, BOND_STROKE_WIDTH, strokeColor);
            createLine(startX - offX, startY - offY, endX - offX, endY - offY, BOND_STROKE_WIDTH, strokeColor);
        } else if (b.type === 'triple') {
            const offX1 = offset * Math.cos(perpAngle);
            const offY1 = offset * Math.sin(perpAngle);
            const offX2 = -offset * Math.cos(perpAngle);
            const offY2 = -offset * Math.sin(perpAngle);
            createLine(startX + offX1, startY + offY1, endX + offX1, endY + offY1, BOND_STROKE_WIDTH, strokeColor);
            createLine(startX, startY, endX, endY, BOND_STROKE_WIDTH, strokeColor);
            createLine(startX + offX2, startY + offY2, endX + offX2, endY + offY2, BOND_STROKE_WIDTH, strokeColor);
        } else if (b.type === 'wedge-up' || b.type === 'wedge-down') {
            createWedge(startX, startY, endX, endY, b.type, strokeColor);
        } else {
            createLine(startX, startY, endX, endY, BOND_STROKE_WIDTH, strokeColor);
        }
    });

    atoms.forEach(a => {
        const shouldDrawSymbol = !(a.elem === 'C' && bondCount[a.id] > 0);

        // 绘制可点击的透明区域（保留所有原子）
        const clickableCircle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        clickableCircle.setAttribute('cx', a.x);
        clickableCircle.setAttribute('cy', a.y);
        clickableCircle.setAttribute('r', ATOM_RADIUS);
        clickableCircle.setAttribute('class', `selectable-atom${selectedAtoms.has(a.id) ? ' selected' : ''}`);
        clickableCircle.setAttribute('data-atom-id', a.id);
        // 选中高亮，使用 HIGHLIGHT_COLOR (淡黄色)
        if (selectedAtoms.has(a.id)) {
            clickableCircle.setAttribute('stroke', 'rgb(255, 255, 0)');
            clickableCircle.setAttribute('stroke-width', 2);
        }
        svg.appendChild(clickableCircle);

        if (shouldDrawSymbol) {
            const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('x', a.x);
            text.setAttribute('y', a.y);
            text.setAttribute('text-anchor', 'middle');
            text.setAttribute('dominant-baseline', 'middle');
            text.setAttribute('fill', colors[a.elem]);
            text.setAttribute('font-size', ATOM_FONT_SIZE);
            text.setAttribute('font-weight', 'bold');
            text.textContent = a.elem;
            text.setAttribute('data-atom-id', a.id);
            svg.appendChild(text);
        }
    });

    if (isSelecting && selectionRect) {
        svg.appendChild(selectionRect);
    }
}

function addAtom(elem, x, y) {
    const newAtom = {
        id: nextAtomId++,
        elem,
        x,
        y
    };
    atoms.push(newAtom);
    draw();
}

function addBond(id1, id2) {
    const existingBondIndex = bonds.findIndex(b =>
        (b.a1 === id1 && b.a2 === id2) || (b.a1 === id2 && b.a2 === id1)
    );

    if (existingBondIndex !== -1) {
        if (currentBond === 'single') {
            const existingBond = bonds[existingBondIndex];
            if (existingBond.type === 'single') {
                existingBond.type = 'double';
            } else if (existingBond.type === 'double') {
                existingBond.type = 'triple';
            } else if (existingBond.type === 'triple') {
                existingBond.type = 'single';
            } else {
                existingBond.type = 'single';
            }
        } else {
            bonds[existingBondIndex].type = currentBond;
        }
    } else {
        bonds.push({
            id: nextBondId++,
            a1: id1,
            a2: id2,
            type: currentBond
        });
    }
    draw();
}

function extendFrom(atomId, mouseX, mouseY) {
    const parentAtom = atoms.find(a => a.id === atomId);
    if (!parentAtom) return;

    const dx = mouseX - parentAtom.x;
    const dy = mouseY - parentAtom.y;
    let angle = Math.atan2(dy, dx);
    if (angle < 0) {
        angle += 2 * Math.PI;
    }

    const bondAngles = [
        Math.PI / 6, Math.PI / 2, 5 * Math.PI / 6,
        7 * Math.PI / 6, 3 * Math.PI / 2, 11 * Math.PI / 6
    ];

    let closestAngle = bondAngles[0];
    let minDiff = Math.abs(angle - closestAngle);
    for (let i = 1; i < bondAngles.length; i++) {
        const diff = Math.abs(angle - bondAngles[i]);
        if (diff < minDiff) {
            minDiff = diff;
            closestAngle = bondAngles[i];
        }
    }

    const dist = 40;
    const newX = parentAtom.x + dist * Math.cos(closestAngle);
    const newY = parentAtom.y + dist * Math.sin(closestAngle);

    const overlappingAtom = atoms.find(a =>
        Math.sqrt(Math.pow(newX - a.x, 2) + Math.pow(newY - a.y, 2)) < ATOM_RADIUS
    );

    if (overlappingAtom) {
        addBond(parentAtom.id, overlappingAtom.id);
    } else {
        const newAtomId = nextAtomId;
        addAtom(currentTool, newX, newY);
        addBond(parentAtom.id, newAtomId);
    }
}

function toggleSelectAtom(atomId) {
    if (selectedAtoms.has(atomId)) {
        selectedAtoms.delete(atomId);
    } else {
        selectedAtoms.add(atomId);
    }
    draw();
}

function toggleSelectBond(bondId) {
    if (selectedBonds.has(bondId)) {
        selectedBonds.delete(bondId);
    } else {
        selectedBonds.add(bondId);
    }
    draw();
}

function clearSelection() {
    selectedAtoms.clear();
    selectedBonds.clear();
    draw();
}

function isInsideSelectionBox(x, y) {
    // 检查点是否在套索框内（稍微放大范围以便于点击）
    const padding = 5;
    return x >= selectionBox.left - padding && x <= selectionBox.right + padding &&
        y >= selectionBox.top - padding && y <= selectionBox.bottom + padding;
}

function getMousePos(e) {
    const rect = svg.getBoundingClientRect();
    return {
        x: e.clientX - rect.left,
        y: e.clientY - rect.top
    };
}

function updateCanvasCursor() {
    svg.classList.remove('move-cursor');
    svg.classList.remove('lasso-cursor');
    if (mode === 'move') {
        svg.classList.add('move-cursor');
        svg.style.cursor = 'move';
    } else if (mode === 'lasso') {
        svg.classList.add('lasso-cursor');
        svg.style.cursor = 'cell';
    } else {
        svg.style.cursor = 'crosshair';
    }
}

/**
 * 删除所有选中的原子和键。
 * 这是处理删除的核心逻辑，供按钮点击和键盘事件调用。
 */
function deleteSelectedElements() {
    const atomsToDelete = new Set(selectedAtoms);
    // bondsToDelete 包含所有被选中的键ID
    const bondsToDelete = new Set(selectedBonds);
    const bondsToRemove = new Set();

    if (atomsToDelete.size === 0 && bondsToDelete.size === 0) return;

    // 1. 如果原子被选中，则与其相连的键也应被删除
    bonds.forEach(b => {
        if (atomsToDelete.has(b.a1) || atomsToDelete.has(b.a2)) {
            bondsToDelete.add(b.id);
        }
    });

    // 2. 过滤掉要删除的键
    bonds = bonds.filter(b => !bondsToDelete.has(b.id));

    // 3. 过滤掉要删除的原子
    atoms = atoms.filter(a => !atomsToDelete.has(a.id));

    // 4. 清除当前选择并重绘
    clearSelection();
}

// --- 鼠标事件监听 ---
svg.addEventListener('mousedown', e => {
    const pos = getMousePos(e);
    dragStartPos = pos;
    const targetElem = e.target;
    const atomId = targetElem.getAttribute('data-atom-id');
    const bondId = targetElem.getAttribute('data-bond-id');

    if (mode === 'lasso') {
        // 检查是否在已有的选中框内点击，且有原子被选中
        if (selectedAtoms.size > 0 && isInsideSelectionBox(pos.x, pos.y)) {
            isDraggingSelection = true;
            // 清除选择框的SVG元素，拖动时不需要显示它
            if (selectionRect && selectionRect.parentNode) {
                selectionRect.parentNode.removeChild(selectionRect);
            }
        } else {
            // 开始新的框选
            isSelecting = true;
            startX = pos.x;
            startY = pos.y;
            // 移除旧的选择框
            if (selectionRect && selectionRect.parentNode) {
                selectionRect.parentNode.removeChild(selectionRect);
            }
            selectionRect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            selectionRect.setAttribute('id', 'selectionRect');
            selectionRect.setAttribute('x', startX);
            selectionRect.setAttribute('y', startY);
            selectionRect.setAttribute('width', 0);
            selectionRect.setAttribute('height', 0);
            svg.appendChild(selectionRect);
            clearSelection(); // 开始新框选时清空旧选择
        }
    } else if (mode === 'move') {
        if (atomId) {
            // 移动模式下点击原子，清除旧选择，并选中当前原子
            clearSelection();
            toggleSelectAtom(parseInt(atomId));
            isDraggingAtom = true;
            dragTargetId = parseInt(atomId);
            const atom = atoms.find(a => a.id === dragTargetId);
            dragStartOffset.x = pos.x - atom.x;
            dragStartOffset.y = pos.y - atom.y;
        } else if (bondId) {
            clearSelection();
            toggleSelectBond(parseInt(bondId));
        } else {
            clearSelection();
        }
    } else { // mode === 'draw'
        if (atomId) {
            isDraggingAtom = true;
            dragTargetId = parseInt(atomId);
        } else if (bondId) {
            const bond = bonds.find(b => b.id === parseInt(bondId));
            if (bond) {
                bond.type = currentBond;
                draw();
            }
        } else if (e.target.tagName === 'svg') {
            addAtom(currentTool, pos.x, pos.y);
            clearSelection();
        }
    }
});

svg.addEventListener('mousemove', e => {
    const pos = getMousePos(e);

    if (isDraggingSelection) {
        const dx = pos.x - dragStartPos.x;
        const dy = pos.y - dragStartPos.y;

        selectedAtoms.forEach(id => {
            const atom = atoms.find(a => a.id === id);
            if (atom) {
                atom.x += dx;
                atom.y += dy;
            }
        });

        // 更新selectionBox
        selectionBox.left += dx;
        selectionBox.right += dx;
        selectionBox.top += dy;
        selectionBox.bottom += dy;

        dragStartPos = pos;
        draw();
        return;
    }

    if (mode === 'lasso' && isSelecting) {
        const width = Math.abs(pos.x - startX);
        const height = Math.abs(pos.y - startY);
        const newX = Math.min(pos.x, startX);
        const newY = Math.min(pos.y, startY);
        selectionRect.setAttribute('x', newX);
        selectionRect.setAttribute('y', newY);
        selectionRect.setAttribute('width', width);
        selectionRect.setAttribute('height', height);
    } else if (mode === 'move' && isDraggingAtom) {
        const atom = atoms.find(a => a.id === dragTargetId);
        if (atom) {
            atom.x = pos.x - dragStartOffset.x;
            atom.y = pos.y - dragStartOffset.y;
            draw();
        }
    } else {
        // 实时悬停效果
        svg.querySelectorAll('.hover-effect').forEach(el => el.remove());
        const targetElem = e.target;
        const atomId = targetElem.getAttribute('data-atom-id');
        const bondId = targetElem.getAttribute('data-bond-id');
        const atomMapForHover = new Map(atoms.map(a => [a.id, a])); // 重新获取 atomMap，防止 draw() 外部调用时不存在

        if (atomId && !selectedAtoms.has(parseInt(atomId))) {
            const atom = atoms.find(a => a.id === parseInt(atomId));
            if (atom) {
                const hoverCircle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                hoverCircle.setAttribute('cx', atom.x);
                hoverCircle.setAttribute('cy', atom.y);
                hoverCircle.setAttribute('r', ATOM_RADIUS);
                hoverCircle.setAttribute('stroke', HOVER_COLOR);
                hoverCircle.setAttribute('stroke-width', 2);
                hoverCircle.setAttribute('fill', 'none');
                hoverCircle.setAttribute('class', 'hover-effect');
                svg.appendChild(hoverCircle);
            }
        } else if (bondId && !selectedBonds.has(parseInt(bondId))) {
            const bond = bonds.find(b => b.id === parseInt(bondId));
            if (bond) {
                const a1 = atomMapForHover.get(bond.a1);
                const a2 = atomMapForHover.get(bond.a2);
                if (!a1 || !a2) return;
                const hoverLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                hoverLine.setAttribute('x1', a1.x);
                hoverLine.setAttribute('y1', a1.y);
                hoverLine.setAttribute('x2', a2.x);
                hoverLine.setAttribute('y2', a2.y);
                hoverLine.setAttribute('stroke', HOVER_COLOR);
                hoverLine.setAttribute('stroke-width', 8);
                hoverLine.setAttribute('class', 'hover-effect');
                svg.appendChild(hoverLine);
            }
        }
    }
});

svg.addEventListener('mouseup', e => {
    const pos = getMousePos(e);
    dragEndPos = pos;
    const targetElem = e.target;
    const atomId = targetElem.getAttribute('data-atom-id');
    const bondId = targetElem.getAttribute('data-bond-id');
    const dragDistance = Math.sqrt(Math.pow(dragEndPos.x - dragStartPos.x, 2) + Math.pow(dragEndPos.y - dragStartPos.y, 2));

    if (isDraggingSelection) {
        isDraggingSelection = false;
        draw();
        return;
    }

    if (mode === 'lasso') {
        isSelecting = false;
        if (selectionRect) {
            svg.removeChild(selectionRect);
            selectionRect = null;
        }

        // 只有在发生了拖动行为时才进行框选，否则视为点击框外重置选择
        if (dragDistance > CLICK_THRESHOLD) {
            const endX = pos.x;
            const endY = pos.y;
            const left = Math.min(startX, endX);
            const top = Math.min(startY, endY);
            const width = Math.abs(startX - endX);
            const height = Math.abs(startY - endY);

            selectionBox.left = left;
            selectionBox.top = top;
            selectionBox.right = left + width;
            selectionBox.bottom = top + height;

            atoms.forEach(a => {
                if (a.x >= left && a.x <= left + width && a.y >= top && a.y <= top + height) {
                    selectedAtoms.add(a.id);
                }
            });
            bonds.forEach(b => {
                const a1 = atoms.find(a => a.id === b.a1);
                const a2 = atoms.find(a => a.id === b.a2);
                if (!a1 || !a2) return;
                if (selectedAtoms.has(a1.id) && selectedAtoms.has(a2.id)) {
                    selectedBonds.add(b.id);
                }
            });
        } else {
            // 如果是点击，且不在选框内，则清空选择
            if (!isInsideSelectionBox(pos.x, pos.y)) {
                clearSelection();
            }
        }
        draw();
    } else if (mode === 'move') {
        isDraggingAtom = false;
        dragTargetId = null;
        // 移动模式下的点击选择已在 mousedown 中处理
        draw();
    } else if (mode === 'draw' && isDraggingAtom) {
        if (dragDistance < CLICK_THRESHOLD) {
            // In draw mode, a click on an atom initiates an extension
        } else {
            const targetAtom = atoms.find(a => a.id !== dragTargetId && Math.sqrt(Math.pow(pos.x - a.x, 2) + Math.pow(pos.y - a.y, 2)) < ATOM_RADIUS);
            if (targetAtom) {
                addBond(dragTargetId, targetAtom.id);
            } else {
                extendFrom(dragTargetId, pos.x, pos.y);
            }
        }
        isDraggingAtom = false;
        dragTargetId = null;
        draw();
    }

    // Re-apply hover effect if needed
    const mouseEvent = new MouseEvent('mousemove', {
        clientX: e.clientX,
        clientY: e.clientY
    });
    svg.dispatchEvent(mouseEvent);
});

svg.addEventListener('mouseleave', () => {
    isSelecting = false;
    isDraggingAtom = false;
    isDraggingSelection = false;
    dragTargetId = null;
    if (selectionRect) {
        svg.removeChild(selectionRect);
        selectionRect = null;
    }
    svg.querySelectorAll('.hover-effect').forEach(el => el.remove());
});

// --- 工具栏按钮逻辑 ---
document.querySelectorAll('#toolbar button').forEach(btn => {
    btn.addEventListener('click', () => {
        const modeButtons = document.querySelectorAll('.mode-btn');

        // 删除按钮有单独的逻辑，不进行模式切换
        if (btn.id === 'deleteBtn') {
            deleteSelectedElements();
            return;
        }

        clearSelection();

        if (btn.classList.contains('mode-btn')) {
            modeButtons.forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            mode = btn.id === 'drawBtn' ? 'draw' : btn.id === 'moveBtn' ? 'move' : 'lasso';
        } else if (btn.classList.contains('atom-btn')) {
            document.querySelectorAll('.atom-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            currentTool = btn.dataset.tool;
            document.getElementById('drawBtn').click();
        } else if (btn.classList.contains('bond-btn')) {
            document.querySelectorAll('.bond-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            currentBond = btn.dataset.bond;
            document.getElementById('drawBtn').click();
        }
        updateCanvasCursor();
    });
});

// 键盘事件监听器，调用统一的删除函数
document.addEventListener('keydown', e => {
    if (e.key === 'Delete' || e.key === 'Backspace') {
        e.preventDefault(); // 阻止浏览器回退等默认行为
        deleteSelectedElements();
    }
});

// 初始设置默认工具
document.getElementById('drawBtn').click();
document.querySelector('.atom-btn[data-tool="C"]').classList.add('active');
document.querySelector('.bond-btn[data-bond="single"]').classList.add('active');
draw();